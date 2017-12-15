% %============================================================================%
% % Duke University                                                            %
% % K. P. Trofatter                                                            %
% % kpt2@duke.edu                                                              %
% %============================================================================%
% Raytrace2d() - finds rays between tx and rx obeying specular reflection.
%
% USAGE:
%   [traces] = Raytrace2d(rx, tx, geometry, depth, debug)
%
% INPUT:
%   [2,t] double  | rx        | [m] receive vertex
%   [2,r] double  | tx        | [m] transmit vertices
%   [1,1] double  | geometry  | scene geometry
%   [2,v] double  | .verts    | [m] vertices
%   [2,e] double  | .edges    | edge vertex indices
%   [1,1] double  | depth     | [#] maximum recursion depth
%   [1,1] logical | debug     | debug plot flag
%
% OUTPUT:
%   [1,n] cell    | traces    | [m] ray traces

function [traces] = Raytrace2d(rx, tx, geometry, depth, debug)
    
DebugPlot(tx, rx, geometry, depth, {});
    % initate traces
    traces = {};
    itrace = 1;
	
    % bake geometry data for faster computation
    geometry = BakeGeometry(geometry);
    
    % trace unoccluded unreflected rx-tx rays
    t = TraceTx(rx, tx, geometry);
    [traces, itrace] = SetTraces(traces, itrace, t);
    
    % list visible geometry
    intervals = GetVisibleEdges(rx, geometry);
    
    % recurse into geometry
    nintervals = numel(intervals);
    for i = 1 : nintervals
        interval = intervals{i};
        if ~isempty(interval)
            t = TraceEdge(rx, tx, geometry, i, interval, depth);
            [traces, itrace] = SetTraces(traces, itrace, t);
        end
    end
    
    % set root of all traces to rx
    for i = 1 : itrace - 1
        traces{i} = [rx, traces{i}];
    end
    
    % prune traces
    traces(itrace : end) = [];
    
    % debug plot
    if debug
        DebugPlot(tx, rx, geometry, depth, traces)
    end
    
end


function [geometry] = BakeGeometry(geometry)
    
    % get geometry
    verts = geometry.verts;
    edges = geometry.edges;
    
    % compute edge line equations
    % x(t) = x0 + u * t
    x0 = verts(:, edges(1, :)); % edge origin
    x1 = verts(:, edges(2, :));
    u = x1 - x0;
    d = sum(u .^ 2) .^ 0.5; % edge length
    u = u ./ repmat(d, [2, 1]); % edge direction
    
    % f(x, y) = A * x + B * y + C = 0
    % f(x, y) gives signed distance to (x, y), (A, B) is normal to line,
    A = x0(2, :) - x1(2, :);
    B = x1(1, :) - x0(1, :);
    C = x0(1, :) .* x1(2, :) - x1(1, :) .* x0(2, :);
    D = (A .^ 2 + B .^ 2) .^ 0.5;
    A = A ./ D;
    B = B ./ D;
    C = C ./ D;
    
    % set geometry
    geometry.x0 = x0;
    geometry.x1 = x1;
    geometry.u = u;
    geometry.d = d;
    geometry.A = A;
    geometry.B = B;
    geometry.C = C;
    
end


function [test] = PointLineIntersect(point, x0, u)
    
    % !!! HARDCODED !!!
    radius = 1.0e-6; % [m] line radius
    
    % compute vectors
    h = point - x0; % hypotenuse vector
    x = dot(h, u); % coordinate of point along edge
    a = x * u;     % adjacent vector
    o = h - a;      % opposite vector
    
    % test if point close enough to line
    test = norm(o) <= radius;
    
end


function [t] = LineLineIntersect(x0, u0, x1, u1)
    
    % parallel test
    d = dot(u0, u1);
    if norm(u0) * norm(u1) - abs(d) < 1.0e-9
        
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


function [a, b] = ShadowInterval(point, x0, x1, x2, x3)
    
    % initialize interval
    a = nan();
    b = nan();
    
    % occluding edge
    u = x1 - x0;
    d = norm(u);
    u = u / d;
    
    % edge end 0
    u0 = x0 - point;
    d0 = norm(u0);
    u0 = u0 ./ d0;
    
    % edge end 1
    u1 = x1 - point;
    d1 = norm(u1);
    u1 = u1 ./ d1;
    
    % edge under test
    ut = x3 - x2;
    dt = norm(ut);
    ut = ut / dt;
    
    % normal pointing into shadow from occlusing edge
    nn = u0 - repmat(dot(u0, u), [2, 1]) .* u;
    nn = nn ./ repmat(sum(nn .^ 2) .^ 0.5, [2,1]);
    
    % test point-occluding line intersect (discards problematic cases)
    if PointLineIntersect(point, x0, u)
        return
    end
    
    % compute intersections
    t  = LineLineIntersect(x0, u , x2, ut);
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
        t1(2) = t(2);
    elseif (-d0 <= t0(1)) && (t0(1) <= 0.0) && (0.0 < t1(1))
        % 0+
        t0(2) = t(2);
    elseif (0.0 < t0(1)) && (t1(1) < -d1)
        % +-
        if dot(nn, ut) >= 0.0
            t1(2) = inf();
        else
            t1(2) = -inf();
        end
    elseif (t0(1) < -d0) && (0.0 < t1(1))
        % -+
        if dot(nn, ut) >= 0.0
            t0(2) = inf();
        else
            t0(2) = -inf();
        end
    elseif (-d0 <= t0(1)) && (t0(1) <= 0.0) && (t1(1) < -d1)
        % 0-
        t0(2) = t(2);
        if dot(nn, ut) >= 0.0
            t1(2) = inf();
        else
            t1(2) = -inf();
        end
    elseif (t0(1) < -d0) && (-d1 <= t1(1)) && (t1(1) <= 0.0)
        % -0
        t0(2) = t(2);
        if dot(nn, ut) >= 0.0
            t1(2) = inf();
        else
            t1(2) = -inf();
        end
    else
        warning('Impossible!?');
    end
    
    % construct interval
    a = t0(2);
    b = t1(2);
    
end


function [intervals] = GetVisibleEdges(point, geometry, imirror, visibility)
    
    % get geometry
    edges = geometry.edges;
    nedges = size(edges, 2);
    x0 = geometry.x0;
    x1 = geometry.x1;
    u = geometry.u;
    d = geometry.d;
    
    % mirror flag
    is_mirror = (nargin() == 4);
    if is_mirror
        x0_mirror = x0(:, imirror);
        u_mirror = u(:, imirror);
    end
    
    % initiate visibility
    intervals = cell(1, nedges);
    for i = 1 : nedges
        intervals{i} = Interval(0.0, d(i));
    end
    
    % edge visibilty loop
    for i = 1 : nedges
        
        % test point-visible line intersect (discards problematic cases)
        if PointLineIntersect(point, x0(:, i), u(:, i))
            intervals{i} = [];
            continue
        end
        
        % get edge interval
        interval = intervals{i};
        
        % occluding edges loop
        for j = 1 : nedges
            
            % test for distinct edges
            if i == j || i == imirror || j == imirror
                continue
            end
            
            % clip occluding edges between image point and mirror
            if is_mirror
                t = LineLineIntersect(x0_mirror, u_mirror, x0(:, j), u(:, j));
                n = u(:, j) - dot(u(:, j), u_mirror) * u_mirror;
                if dot(n, u(:, j)) >= 0.0
                    b = -inf();
                else
                    b = inf();
                end
                
            end
            
            % test for empty interval
            if isempty(interval)
                break
            end
            
            % compute shadow interval
            [a, b] = ShadowInterval(point, x0(:, j), x1(:, j), x0(:, i), x1(:, i));
            
            % update visible interval
            if ~isnan(a) && ~isnan(b)
                interval = IntervalDifference(interval, a, b);
            end
            
        end
        
        % set inteval
        intervals{i} = interval;
        
    end
    
end


function [traces] = TraceTx(rx, tx, geometry, imirror, interval)
    
    % mirror flag
    is_mirror = (nargin() == 5);
    
    % get mirror edge
    if is_mirror
        mirror_x0 = geometry.x0(:, imirror);
        mirror_u = geometry.u(:, imirror);
    end
    
    % count things
    ntx = size(tx, 2);
    nedges = size(geometry.edges, 2);
    
    % construct rx-to-tx vector line equations
    % x(t) = x0 + u * t
    x0 = repmat(rx, [1, ntx]);
    u = tx - x0;
    d = sum(u .^ 2) .^ 0.5;
    u = u ./ repmat(d, [2, 1]);
    
    % initiate ray traces
    traces = cell(1, ntx);
    itrace = 1;
    
    % for all tx
    for i = 1 : ntx
        
        % set occlusion minimum along rx-tx line, test mirror interval mask
        if is_mirror
            if PointLineIntersect(tx(:, i), mirror_x0, mirror_u)
                continue
            end
            % test for visibility through mask
            t = LineLineIntersect(mirror_x0, mirror_u, x0(:, i), u(:, i));
            is_visible = IntervalElementOf(interval, t(1)) ...
                && (0.0 <= t(2)) && (t(2) <= d(i));
            
            % not visible, skip to next tx
            if ~is_visible()
                continue
            end
            
            % set minimum at mirror surface (edges behind mirror do not occlude)
            % NOTE : this is works in conjunction with the edge test
            minimum = t(2);
            
        else
            
            % set visibility flag
            is_visible = true();
            
            % set occlusion minimum along rx-tx line at rx point
            minimum = 0.0;
            
        end
        
        % test for occluding edges
        for j = 1 : nedges
            
            % skip mirror edge (we just took care of it)
            if is_mirror && (j == imirror)
                continue
            end
            
            % toss out problematic cases
            edge_x0 = geometry.x0(:, j);
            edge_u = geometry.u(:, j);
            if PointLineIntersect(tx(:, i), edge_x0, edge_u)
                continue
            end
            
            % find intersection between rx-tx line and edge line
            t = LineLineIntersect(edge_x0, edge_u, x0(:, i), u(:, i));
            
            % test if intersection on edge, and edge is occluding
            edge_d = geometry.d(j);
            if (0.0 <= t(1)) && (t(1) <= edge_d) ...
            && (minimum <= t(2)) && (t(2) <= d(i))
                
                % occluded!
                is_visible = false();
                break
                
            end
            
        end
        
        % set trace
        if is_visible
            traces{itrace} = tx(:, i);
            itrace = itrace + 1;
        end
        
    end
    
    % prune traces
    traces(itrace : end) = [];

end


function [traces] = TraceEdge(point, tx, geometry, imirror, interval, depth)
    
    % initiate traces
    traces = {};
    itrace = 1;
    
    % mirror point
    u_mirror = geometry.u(:, imirror);
    x0_mirror = geometry.x0(:, imirror);
    r = point - x0_mirror;
    image = 2.0 * dot(r, u_mirror) * u_mirror - r + x0_mirror;
    
    % trace tx
    t = TraceTx(image, tx, geometry, imirror, interval);
    [traces, itrace] = SetTraces(traces, itrace, t);
    
    % recurse
    if depth > 0
        
        % get edge visibility
        intervals = GetVisibleEdges(image, geometry, imirror, interval);
        
        % trace edges
        nintervals = numel(intervals);
        for i = 1 : nintervals
            t = TraceEdge(point, tx, geometry, i, intervals{i}, depth - 1);
            [traces, itrace] = SetTraces(traces, itrace, t);
        end
        
    end
    
    % push reflected point to traces
    for i = 1 : itrace - 1
        x0 = traces{i}(:, 1);
        u = image - x0;
        u = u / norm(u);
        t = LineLineIntersect(x0_mirror, u_mirror, x0, u);
        reflect = x0 + u * t(2);
        traces{i} = [reflect, traces{i}];
    end
    
    % prune traces
    traces(itrace : end) = [];
    
end


function [traces, itrace] = SetTraces(traces, itrace, t)
    
    % allocate
    if itrace > numel(traces)
        traces = [traces, cell(1, 64)];
    end
    
    % set
    n = numel(t);
    traces(itrace + (0 : n - 1)) = t;
    itrace = itrace + n;
    
end


function [interval] = Interval(a, b)
    interval = [a; b];
end


function [test] = IntervalElementOf(interval, a)
    
    nintervals = size(interval, 2);
    for i = 1 : nintervals
        if interval(1, i) <= a && a <= interval(2, i)
            test = true();
            return
        end
    end
    test = false();
    
end


function [intersect] = IntervalIntersect(interval, a, b)
    
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
        
        % compute set difference
        if (b < c) || (d < a)
            % ab cd, cd ab
        elseif (a <= c) && (d <= b)
            % a cd b
            intersect(:, idiff) = [c; d];
            idiff = idiff + 1;
        elseif (c <= a) && (b <= d)
            % c ab d
            intersect(:, idiff) = [a; b];
            idiff = idiff + 1;
        elseif (a <= c) && (b <= d)
            % a c b d
            intersect(:, idiff) = [c; b];
            idiff = idiff + 1;
        elseif (c <= a) && (d <= b)
            % c a d b
            intersect(:, idiff) = [a; d];
            idiff = idiff + 1;
        end
            
    end
    
    % prune
    intersect(:, idiff : end) = [];
    
end


function [diff] = IntervalDifference(interval, a, b)
    
    % order points
    if a > b
        c = a;
        a = b;
        b = c;
    end
    
    % allocate new interval
    nintervals = size(interval, 2);
    diff = zeros(2, 2 * nintervals);
    idiff = 1;
    
    for i = 1 : nintervals
        
        % get subinterval
        c = interval(1, i);
        d = interval(2, i);
        
        % compute set difference
        if (b <= c) || (d <= a)
            % ab cd, cd ab
            diff(:, idiff) = [c; d];
            idiff = idiff + 1;
        elseif (a <= c) && (d <= b)
            % a cd b
        elseif (c < a) && (b < d)
            % c ab d
            diff(:, idiff + 0) = [c; a];
            diff(:, idiff + 1) = [b; d];
            idiff = idiff + 2;
        elseif (a <= c) && (b <= d)
            % a c b d
            diff(:, idiff) = [b; d];
            idiff = idiff + 1;
        elseif (c <= a) && (d <= b)
            % c a d b
            diff(:, idiff) = [c; a];
            idiff = idiff + 1;
        end
            
    end
    
    % prune
    diff(:, idiff : end) = [];
    
end


function [] = DebugPlot(tx, rx, geometry, depth, traces)
    
    % figure
    fh = figure(1);
    clf(fh);
    ah = axes('Parent', fh);
    
    % draw begin
    hold(ah, 'on');
    
    % tx, rx
    scatter(ah, tx(1, :), tx(2, :), 'r', 'filled');
    scatter(ah, rx(1, :), rx(2, :), 'b', 'filled');
    
    % verts
    verts = geometry.verts;
    scatter(ah, verts(1, :), verts(2, :), 'k', 'filled');
    
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
        'Color', 'k', 'LineWidth', 2.0, 'Parent', ah);
    
    % rays
    nrays = numel(traces);
    for i = 1 : nrays
        line('XData', traces{i}(1, :), 'YData', traces{i}(2, :), 'Color', 'r', 'Parent', ah);
    end
    
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


%==============================================================================%
%                                                                              %
%                                                                              %
%                                                                              %
%==============================================================================%
