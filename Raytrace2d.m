% %============================================================================%
% % Duke University                                                            %
% % K. P. Trofatter                                                            %
% % kpt2@duke.edu                                                              %
% %============================================================================%
% Raytrace2d() - finds rays between tx and rx obeying specular reflection.
%
% USAGE:
%   [traces] = Raytrace2d(rx, tx, geometry, depth, debug?)
%
% INPUT:
%   [2,r] double  | rx        | [m] receive vertex
%   [2,t] double  | tx        | [m] transmit vertices
%   [1,1] double  | geometry  | scene geometry
%   [2,v] double  | .verts    | [m] vertices
%   [2,e] double  | .edges    | edge vertex indices
%   [1,1] double  | depth     | [#] maximum recursion depth
%   [1,1] char    | debug     | debug plot {'solid', 'bounce', 'distance'}
%
% OUTPUT:
%   [1,n] cell    | traces    | [m] ray traces ([2,?] double [rx, ... ,tx])
%
% NOTES:
%   + currently only supports edge geometry

function [traces] = Raytrace2d(rx, tx, geometry, depth, debug)
    
    % default debug type
    if ~exist('debug', 'var')
        debug = [];
    end
    
    % debug plot
    if debug
        
        % figure and axes
        fh = DarkFigure(1);
        clf(fh);
        ah = DarkAxes(fh);
        
        % plot scene
        DebugPlotScene(ah, tx, rx, geometry, depth);
        
    end
	
    % bake reusable geometry computations
    geometry = BakeGeometry(geometry);
    
    % recursively trace
    traces = TraceInto(rx, tx, geometry, [], [], depth);
    
    % set root of all traces to rx
    for i = 1 : numel(traces)
        traces{i} = [rx, traces{i}];
    end
    
    % plot traces
    if debug
        DebugPlotTraces(ah, traces, depth, debug);
    end
    
end


function [geometry] = BakeGeometry(geometry)
    
    % get geometry
    verts = geometry.verts;
    edges = geometry.edges;
    
    % edge vector line equations
    % x(t) = x0 + u * t
    x0 = verts(:, edges(1, :)); % edge begin
    x1 = verts(:, edges(2, :)); % edge end
    u = x1 - x0;
    d = sqrt(sum(u .^ 2));      % edge distance
    u = u ./ repmat(d, [2, 1]); % edge unit direction
    
    % set geometry
    geometry.x0 = x0;
    geometry.x1 = x1;
    geometry.u = u;
    geometry.d = d;
    
end


function [traces] = TraceInto(rx, tx, geometry, imirror, visibility, depth)
    
    % initiate traces
    [traces, itrace] = Traces();
    
    % trace terminal tx
    t = TraceTx(rx, tx, geometry, imirror, visibility);
    [traces, itrace] = TracesSet(traces, itrace, t);
    
    % recursively trace geometry
    if depth > 0
        
        % list visible edges
        intervals = GetVisibleEdges(rx, geometry, imirror, visibility);
        
        % trace visible edges
        for i = 1 : numel(intervals)
            if ~isempty(intervals{i})
                t = TraceEdge(rx, tx, geometry, i, intervals{i}, depth);
                [traces, itrace] = TracesSet(traces, itrace, t);
            end
        end
        
        % add other geometry handlers here, i.e. circles, parabolas, etc.
        
    end
    
    % prune unused traces
    traces(itrace : end) = [];
    
end


% tidy up
function [traces] = TraceTx(rx, tx, geometry, imirror, visibility)
    
    % mirror flag
    is_mirror = ~isempty(imirror);
    
    % get mirror edge
    if is_mirror
        mirror_x0 = geometry.x0(:, imirror);
        mirror_u = geometry.u(:, imirror);
        mirror_d = geometry.d(imirror);
    end
    
    % count things
    ntx = size(tx, 2);
    nedges = size(geometry.edges, 2);
    
    % construct rx->tx vector line equations
    % x(t) = x0 + u * t
    x0 = repmat(rx, [1, ntx]);
    u = tx - x0;
    d = sqrt(sum(u .^ 2));
    u = u ./ repmat(d, [2, 1]);
    
    % initiate traces
    [traces, itrace] = Traces();
    
    % for all tx
    for i = 1 : ntx
        
        % set minimum occlusion distance along rx->tx line, initiate visibility
        if is_mirror
            
            % tx embedded in mirrors should not reflect from them
            if PointLineIntersect(tx(:, i), mirror_x0, mirror_u)
                continue
            end
            
            % test visibility through mask
            t = LineLineIntersect(mirror_x0, mirror_u, x0(:, i), u(:, i));
            n = [mirror_u(2); -mirror_u(1)];
            is_visible = (0.0 <= t(1)) && (t(1) <= mirror_d) ...
                && IntervalElementOf(visibility, t(1)) ...
                && dot(n, tx(:, i) - mirror_x0) * dot(n, rx - mirror_x0) < 0.0;
            
            % not visible
            if ~is_visible()
                continue
            end
            
            % set minimum at mirror surface (walls behind mirror do not occlude)
            minimum = t(2);
            
        else
            
            % set visibility flag
            is_visible = true();
            
            % set occlusion minimum along rx->tx line at rx point
            minimum = 0.0;
            
        end
        
        % for all wall edges
        for j = 1 : nedges
            
            % skip mirror edge (we just took care visibility)
            if is_mirror && (j == imirror)
                continue
            end
            
            % get wall geometry
            wall_x0 = geometry.x0(:, j);
            wall_u = geometry.u(:, j);
            wall_d = geometry.d(j);
            
            % tx embedded in walls should not be occluded by them
            if PointEdgeIntersect(tx(:, i), wall_x0, wall_u, wall_d)
                continue
            end
            
            % find intersection between wall line and rx->tx line
            t = LineLineIntersect(wall_x0, wall_u, x0(:, i), u(:, i));
            
            % test if rx->tx ray intersects wall
            if (0.0 <= t(1)) && (t(1) <= wall_d) ...
            && (minimum <= t(2)) && (t(2) <= d(i))
                
                % occluded!
                is_visible = false();
                break
                
            end
            
        end
        
        % set trace
        if is_visible
            [traces, itrace] = TracesSet(traces, itrace, {tx(:, i)});
        end
        
    end
    
    % prune traces
    traces(itrace : end) = [];

end


function [traces] = TraceEdge(point, tx, geometry, imirror, visibility, depth)
    
    % compute image point
    x0_mirror = geometry.x0(:, imirror);
    u_mirror = geometry.u(:, imirror);
    r = point - x0_mirror;
    image = x0_mirror + 2.0 * dot(r, u_mirror) * u_mirror - r;
    
    % recursively trace tx
    traces = TraceInto(image, tx, geometry, imirror, visibility, depth - 1);
    
    % push reflect point to traces
    for i = 1 : numel(traces)
        
        % compute reflect point
        x0 = traces{i}(:, 1);
        u = image - x0;
        u = u / norm(u);
        t = LineLineIntersect(x0, u, x0_mirror, u_mirror);
        reflect = x0 + u * t(1);
        
        % push
        traces{i} = [reflect, traces{i}];
        
    end
    
end


function [traces, itrace] = Traces()
    traces = {};
    itrace = 1;
end


function [traces, itrace] = TracesSet(traces, itrace, t)
    
    % !!! HARDCODED !!!
    guess = 64; % [#] number of traces to allocate
    
    % count number of traces to set
    n = numel(t);
    
    % allocate
    if itrace + n > numel(traces)
        traces = [traces, cell(1, max(n, guess))];
    end
    
    % set
    traces(itrace + (0 : n - 1)) = t;
    itrace = itrace + n;
    
end


% solid logic, but tidy up a bit
function [intervals] = GetVisibleEdges(rx, geometry, imirror, visibility)
    
    % mirror flag
    is_mirror = ~isempty(imirror);
    
    % get geometry
    nedges = size(geometry.edges, 2);
    x0 = geometry.x0;
    x1 = geometry.x1;
    u = geometry.u;
    d = geometry.d;
    if is_mirror
        x0_mirror = x0(:, imirror);
        u_mirror = u(:, imirror);
    end
    
    % initiate visibility
    intervals = cell(1, nedges);
    
    % for all test edges
    for i = 1 : nedges
        
        % a mirror cannot directly see itself
        % a point embedded in an edge does not see the edge
        if is_mirror && (i == imirror) ...
        || PointEdgeIntersect(rx, x0(:, i), u(:, i), d(i))
            continue
        end
        
        % initiate edge visibility interval
        if is_mirror
            
            % initiate interval with regions visible through mirror
            nsubintervals = size(visibility, 2);
            interval = zeros(2, nsubintervals);
            for j = 1 : nsubintervals
                x2 = x0_mirror + u_mirror * visibility(1, j);
                x3 = x0_mirror + u_mirror * visibility(2, j);
                [a, b] = EdgeShadowInterval(rx, x2, x3, x0(:, i), x1(:, i));
                interval(1 : 2, j) = [a; b];
            end
            interval(:, isnan(interval(1, :))) = [];
            if ~isempty(interval)
                interval = IntervalIntersect(interval, 0.0, d(i));
            end
            
            % compute direction from mirror to shadow
            n_mirror = x0_mirror - rx;
            n_mirror = n_mirror - dot(n_mirror, u_mirror) * u_mirror;
            
        else
            
            % initiatlize to length of edge
            interval = Interval(0.0, d(i));
            
        end
        
        % for all wall edges
        for j = 1 : nedges
            
            % test for self-occlusion, mirror
            if (i == j) || is_mirror && (j == imirror)
                continue
            end
            
            % test for empty interval
            if isempty(interval)
                break
            end
            
            % rx embedded in edges are not occluded by those edges
            if PointEdgeIntersect(rx, x0(:, j), u(:, j), d(j))
                continue
            end
            
            % clip walls on image side of mirror line
            if is_mirror
                a = 0.0;
                b = d(j);
                c1 = (dot(x0(:, j) - x0_mirror, n_mirror) < 0.0);
                c2 = (dot(x1(:, j) - x0_mirror, n_mirror) < 0.0);
                if c1 && c2
                    % fully clipped
                    continue
                elseif (c1 && ~c2) || (~c1 && c2)
                    % clipped accross mirror
                    t = LineLineIntersect(x0(:, j), u(:, j), x0_mirror, u_mirror);
                    if c1
                        a = t(1);
                    else
                        b = t(1);
                    end
                end
                cut_x0 = x0(:, j) + u(:, j) * a;
                cut_x1 = x0(:, j) + u(:, j) * b;
            else
                % default points
                cut_x0 = x0(:, j);
                cut_x1 = x1(:, j);
            end
            
            % compute wall shadow interval on test line
            [a, b] = EdgeShadowInterval(rx, cut_x0, cut_x1, x0(:, i), x1(:, i));
            
            % update visible interval on test line
            if ~isnan(a) && ~isnan(b)
                interval = IntervalDifference(interval, a, b);
            end
            
        end
        
        % set inteval
        intervals{i} = interval;
        
    end
    
end


function [test] = PointLineIntersect(point, x0, u)
    
    % !!! HARDCODED !!!
    radius = 1.0e-6; % [m] line radius
    
    % compute vectors
    h = point - x0; % hypotenuse vector
    x = dot(h, u);  % coordinate of point along edge
    a = x * u;      % adjacent vector
    o = h - a;      % opposite vector
    
    % test if point close enough to line
    test = norm(o) <= radius;
    
end


function [test] = PointEdgeIntersect(point, x0, u, d)
    x = dot(point, u);
    test = PointLineIntersect(point, x0, u) && (0.0 <= x) && (x <= d);
end


function [t] = LineLineIntersect(x0, u0, x1, u1)
    
    % !!! HARDCODED !!!
    epsilon = 1.0e-8; % [#] parallel test threshold
    
    % parallel test
    d = dot(u0, u1);
    if norm(u0) * norm(u1) - abs(d) < epsilon
        
        % parallel, intersection at infinity (even when overlapping)
        t = inf(2, 1);
        if d < 0.0
            t(2) = -inf();
        end
        
    else
        
        % compute intersection along each line
        t = [u0, -u1] \ (x1 - x0);
        
    end
    
end


% works, but could handle corner cases better
% try implementing the implict strategy and benchmark both
function [a, b] = EdgeShadowInterval(point, x0, x1, x2, x3)
    
    % initialize interval
    a = nan();
    b = nan();
    
    if x0 == x1
        return
    end
    
    % wall edge vector line equation
    uw = x1 - x0;
    dw = norm(uw);
    uw = uw / dw;
    
    % edge end 0
    u0 = x0 - point;
    d0 = norm(u0);
    u0 = u0 ./ d0;
    
    % edge end 1
    u1 = x1 - point;
    d1 = norm(u1);
    u1 = u1 ./ d1;
    
    % test edge
    ut = x3 - x2;
    dt = norm(ut);
    ut = ut / dt;
    
    % normal pointing from wall into shadow
    ws = u0 - repmat(dot(u0, uw), [2, 1]) .* uw;
    ws = ws ./ repmat(sum(ws .^ 2) .^ 0.5, [2,1]);
    
    % test point-wall line intersect (discards problematic cases)
    if PointLineIntersect(point, x0, uw)
        return
    end
    
    % compute intersections
    tw = LineLineIntersect(x0, uw, x2, ut);
    t0 = LineLineIntersect(x0, u0, x2, ut);
    t1 = LineLineIntersect(x1, u1, x2, ut);
    
    % build occluding interval
    if (0.0 < t0(1)) && (0.0 < t1(1))
        % ++ no tweaks
    elseif (-d0 <= t0(1)) && (t0(1) <= 0.0) ...
            && (-d1 <= t1(1)) && (t1(1) <= 0.0)
        % 00
        return
    elseif (t0(1) < -d0) && (t1(1) < -d1)
        % --
        return
    elseif (0.0 < t0(1)) && (-d1 <= t1(1)) && (t1(1) <= 0.0)
        % +0
        t1(2) = tw(2);
    elseif (-d0 <= t0(1)) && (t0(1) <= 0.0) && (0.0 < t1(1))
        % 0+
        t0(2) = tw(2);
    elseif (0.0 < t0(1)) && (t1(1) < -d1)
        % +-
        if dot(ws, ut) >= 0.0
            t1(2) = inf();
        else
            t1(2) = -inf();
        end
    elseif (t0(1) < -d0) && (0.0 < t1(1))
        % -+
        if dot(ws, ut) >= 0.0
            t0(2) = inf();
        else
            t0(2) = -inf();
        end
    elseif (-d0 <= t0(1)) && (t0(1) <= 0.0) && (t1(1) < -d1)
        % 0-
        t0(2) = tw(2);
        if dot(ws, ut) >= 0.0
            t1(2) = inf();
        else
            t1(2) = -inf();
        end
    elseif (t0(1) < -d0) && (-d1 <= t1(1)) && (t1(1) <= 0.0)
        % -0
        t0(2) = tw(2);
        if dot(ws, ut) >= 0.0
            t1(2) = inf();
        else
            t1(2) = -inf();
        end
    else
        warning('Impossible!?');
    end
    
    % construct interval
    if t0(2) <= t1(2)
        a = t0(2);
        b = t1(2);
    else
        a = t1(2);
        b = t0(2);
    end
    
end


function [interval] = Interval(a, b)
    interval = [a; b]; % [smaller; larger];
end


function [test] = IntervalElementOf(interval, a)
    
    % test if element of subintervals
    for i = 1 : size(interval, 2)
        if interval(1, i) <= a && a <= interval(2, i)
            test = true();
            return
        end
    end
    test = false();
    
end


function [intersect] = IntervalIntersect(interval, a, b)
    
    % point intervals leave nothing behind
    if a == b
        intersect = zeros(2, 0);
        return
    end
    
    % order points
    if a > b
        c = a;
        a = b;
        b = c;
    end
    
    % allocate new interval
    nintervals = size(interval, 2);
    intersect = zeros(2, nintervals);
    idiff = 1;
    
    for i = 1 : nintervals
        
        % get subinterval
        c = interval(1, i);
        d = interval(2, i);
        
        % compute set intersection
        if (b <= c) || (d <= a)
            % ab c=d, c=d ab
        elseif (a <= c) && (d <= b)
            % a c=d b
            intersect(:, idiff) = [c; d];
            idiff = idiff + 1;
        elseif (c <= a) && (b <= d)
            % c=ab=d
            intersect(:, idiff) = [a; b];
            idiff = idiff + 1;
        elseif (a <= c) && (b <= d)
            % a c=b=d
            intersect(:, idiff) = [c; b];
            idiff = idiff + 1;
        elseif (c <= a) && (d <= b)
            % c=a=d b
            intersect(:, idiff) = [a; d];
            idiff = idiff + 1;
        end
        
    end
    
    % prune
    intersect(:, idiff : end) = [];
    
end


function [diff] = IntervalDifference(interval, a, b)
    
    % point intervals take nothing away
    if a == b
        diff = interval;
        return
    end
    
    % order points
    if a > b
        c = a;
        a = b;
        b = c;
    end
    
    % allocate new interval
    nintervals = size(interval, 2);
    diff = zeros(2, nintervals + 1);
    idiff = 1;
    
    for i = 1 : nintervals
        
        % get subinterval
        c = interval(1, i);
        d = interval(2, i);
        
        % compute set difference
        if (b <= c) || (d <= a)
            % ab c=d, c=d ab
            diff(:, idiff) = [c; d];
            idiff = idiff + 1;
        elseif (a <= c) && (d <= b)
            % a c=d b
        elseif (c < a) && (b < d)
            % c=ab=d
            diff(:, idiff + 0) = [c; a];
            diff(:, idiff + 1) = [b; d];
            idiff = idiff + 2;
        elseif (a <= c) && (b <= d)
            % a c=b=d
            diff(:, idiff) = [b; d];
            idiff = idiff + 1;
        elseif (c <= a) && (d <= b)
            % c=a=d b
            diff(:, idiff) = [c; a];
            idiff = idiff + 1;
        end
        
    end
    
    % prune
    diff(:, idiff : end) = [];
    
end


function [ah] = DebugPlotScene(ah, tx, rx, geometry, depth)
    
    % draw begin
    hold(ah, 'on');
    
    % tx, rx
    scatter(ah, tx(1, :), tx(2, :), 'r', 'filled');
    scatter(ah, rx(1, :), rx(2, :), 'b', 'filled');
    
    % verts
    verts = geometry.verts;
%    scatter(ah, verts(1, :), verts(2, :), 'k', 'filled');
    
    % edges
    edges = geometry.edges;
    nedges = size(edges, 2);
    e = zeros(2, nedges * 3);
    for i = 1 : nedges
        e1 = verts(:, edges(1, i));
        e2 = verts(:, edges(2, i));
        indices = 3 * i + (-2 : 0);
        e(:, indices) = [e1, e2, nan(2, 1)];
    end
    line('Xdata', e(1, :), 'Ydata', e(2, :), ...
        'Color', [0.6, 0.6, 0.6], 'LineWidth', 2.0, 'Parent', ah);
    
    % draw end
    hold(ah, 'off');
    
    % format
    axis(ah, 'equal');
    grid(ah, 'on');
    
    % notate
    xlabel(ah, 'x[m]');
    ylabel(ah, 'y[m]');
    str = sprintf('Raytrace Depth = %i', depth);
    title(ah, str);
    
end


function [] = DebugPlotTraces(ah, traces, depth, debug)
    
    % draw begin
    hold(ah, 'on');
    
    % traces
    ntraces = numel(traces);
    for i = 1 : ntraces
        switch debug
        case 'solid'
            % solid colored lines
            line('XData', traces{i}(1, :), 'YData', traces{i}(2, :), ...
                'Color', 'g', 'Parent', ah);
        
        case 'bounce'
            % colormap
            colormap(flip(parula(depth + 1)));
            ch = colorbar(ah, 'Color', 'w');
            title(ch, 'ray[#]', 'Color', 'w');
            caxis([1, depth + 1]);
            trace = traces{i};
            % kludge to draw colored lines
            x = trace(1, :);
            y = trace(2, :);
            z = zeros(size(x));
            c = flip((1 : size(trace, 2)) - 1, 2);
            surface([x; x], [y; y], [z; z], [c; c],...
                'FaceColor', 'no',...
                'EdgeColor', 'flat',...
                'LineWidth', 1.0);
        
        case 'distance'
            % colormap
            colormap(flip(parula()));
            ch = colorbar(ah, 'Color', 'w');
            title(ch, 'ray[m]', 'Color', 'w');
            trace = traces{i};
            % kludge to draw colored lines
            x = trace(1, :);
            y = trace(2, :);
            z = zeros(size(x));
            c = trace(:, 1 : end - 1) - trace(:, 2 : end);
            c = sqrt(sum(c .^ 2));
            c = [cumsum(c, 'reverse'), 0.0]; % yum
            surface([x; x], [y; y], [z; z], [c; c],...
                'FaceColor', 'no',...
                'EdgeColor', 'interp',...
                'LineWidth', 1.0);
            
        end
        
    end
    
    % draw end
    hold(ah, 'off');
    
end


%==============================================================================%
%                                                                              %
%                                                                              %
%                                                                              %
%==============================================================================%
