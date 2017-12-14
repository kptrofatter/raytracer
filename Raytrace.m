% %============================================================================%
% % Duke University                                                            %
% % K. P. Trofatter                                                            %
% % kpt2@duke.edu                                                              %
% %============================================================================%
% Raytrace() - finds rays between tx and rx obeying specular reflection.
%
% USAGE:
%   [rays] = Raytrace(tx, rx, geometry, depth, debug)
%
% INPUT:
%   [2,t] double  | tx        | [m] transmit vertices
%   [2,r] double  | rx        | [m] receive vertex
%   [1,1] double  | geometry  | scene geometry
%   [2,v] double  | .verts    | [m] vertices
%   [2,e] double  | .edges    | edge vertex indices
%   [1,1] double  | depth     | [#] maximum recursion depth
%   [1,1] logical | debug     | debug plot flag
%
% OUTPUT:
%   [1,n] cell    | rays      | [m] ray traces

    % rays = raytrace(tx, rx, geometry, depth)
    %   rays = []
    %   for all visible tx (w.r.t. rx)
    %     append to rays (trivial case)
    %     record tx
    %   for all mirror = visible geometry (w.r.t. rx)
    %     append to rays (recursive start)
    %     record trace(tx, rx, mirror, geometry, depth)
    %   for all rays
    %     record rx
    
    % rays = trace(tx, point, mirror, geometry, depth)
    %   rays = []
    %   image = mirrored point
    %   for all visible tx (w.r.t. image thru mirror)
    %     append to rays (base case)
    %     record tx
    %   if search depth not exceeded
    %     for all mirror2 = visible geometry (w.r.t. image thru mirror)
    %       append to rays (recursive case)
    %       record trace(tx, image, mirror2, geometry, depth - 1)
    %   for all rays
    %     record reflection point
    

% point to source
% point to edges => solves initiate visiblity problem
% image point to source
% image point to edges

% point-line intersection
% line-line intersection


function [rays] = Raytrace(tx, rx, geometry, depth, debug)
    
    rays = [];
    
    % augment geometry
    [geometry] = AugmentGeometry(geometry);
    
    DebugPlot(tx, rx, geometry, depth, rays);

    [visibility] = Visibility(rx, geometry);
    
%     if debug
%         DebugPlot(tx, rx, geometry, depth, rays)
%     end
    
end


function [t] = LineLineIntersect(v0, n0, v1, n1)
    
    % parallel test
    d = dot(n0, n1);
    if abs(abs(d) - 1.0) < 1.0e-8
        
        % parallel, intersection at infinity (even when overlapping)
        t = inf(2, 1);
        if d < 0.0
            t(2) = -inf();
        end
        
    else
        
        % compute intersection along each line
        t = [n0, -n1] \ (v1 - v0);
        
    end
    
end


function [a, b] = LineHalfPlaneIntersect(v0, n0, v1, n1, np)
    
    % get line-line intersection
    t = LineLineIntersect(v0, n0, v1, n1);
    
    % build interval
    a = t(1);
    if dot(n0, n2) >= 0.0
        b = inf();
    else
        b = -inf();
    end
    
end


function [a, b] = EdgeShadowInterval(point, u0, u1, v0, v1)
    
    % compute bounding line equations
    
    
    
    % x(t) = t * n + v0
    n = v1 - v0;
    d = sum(n .^ 2) .^ 0.5; % edge length
    n = n ./ repmat(d, [2, 1]); % edge direction
    
    % compute line equations
    A = v0(2, :) - v1(2, :);
    B = v1(1, :) - v0(1, :);
    C = v0(1, :) .* v1(2, :) - v1(1, :) .* v0(2, :);
    D = (A .^ 2 + B .^ 2) ^ 0.5;
    A = A ./ D;
    B = B ./ D;
    C = C ./ D;
    
    

end


function [geometry] = AugmentGeometry(geometry)
    
    % get geometry
    verts = geometry.verts;
    edges = geometry.edges;
    
    % compute edge line equations
    % x(t) = t * n + v0
    v0 = verts(:, edges(1, :)); % edge origin
    v1 = verts(:, edges(2, :));
    n = v1 - v0;
    d = sum(n .^ 2) .^ 0.5; % edge length
    n = n ./ repmat(d, [2, 1]); % edge direction
    
    % f(x, y) = A * x + B * y + C = 0
    % f(x, y) gives signed distance to (x, y), (A, B) is normal to line,
    A = v0(2, :) - v1(2, :);
    B = v1(1, :) - v0(1, :);
    C = v0(1, :) .* v1(2, :) - v1(1, :) .* v0(2, :);
    D = (A .^ 2 + B .^ 2) ^ 0.5;
    A = A ./ D;
    B = B ./ D;
    C = C ./ D;
    
    % set geometry
    geometry.v0 = v0;
    geometry.v1 = v1;
    geometry.n = n;
    geometry.d = d;
    geometry.A = A;
    geometry.B = B;
    geometry.C = C;
    
end

function [test] = PointLineIntersect(point, v0, n0)
    
    % !!! HARDCODED !!!
    radius = 1.0e-5; % [m] line radius
    
    % compute vectors
    h = point - v0; % hypotenuse vector
    x = dot(h, n0); % coordinate of point along edge
    a = x * n0;     % adjacent vector
    o = h - a;      % opposite vector
    
    % test if point close enough to line
    test = norm(o) <= radius;
    
end




function [visibility] = Visibility(point, geometry)

    
    
    
    % compute point-to-edge-end line equations
    P = repmat(point, [1, nverts]);
    % edge end 0
    n0 = v0 - P;
    d0 = sum(n0 .^ 2) .^ 0.5;
    n0 = n0 ./ repmat(d0, [2, 1]);
    % edge end 1
    n1 = v1 - P;
    d1 = sum(n1 .^ 2) .^ 0.5;
    n1 = n1 ./ repmat(d1, [2, 1]);
    
    nn = n0 - repmat(dot(n0, n), [2, 1]) .* n;
    nn = nn ./ repmat(sum(nn .^ 2) .^ 0.5, [2,1]);
    
    % initiate visibility
    visibility = cell(1, nedges);
    for i = 1 : nedges
        visibility{i} = Interval(0.0, d(i));
    end
    
    % edge visibilty loop
    for i = 1 : nedges
        
        % test point-visible line intersect (discards problematic cases)
        if PointLineIntersect(point, v0(:, i), n(:, i))
            visibility{i} = [];
            continue
        end
        
        % get edge interval
        interval = visibility{i};
        
        % occluding edges loop
        for j = 1 : nedges
            
            % test for self occlusion
            if i == j
                continue
            end
            
            % test for empty interval
            if isempty(interval)
                break
            end
            
            % test point-occluding line intersect (discards problematic cases)
            if PointLineIntersect(point, v0(:, j), n(:, j))
                continue
            end
            
            % compute intersections
            t  = LineLineIntersect(v0(:, j), n(:, j) , v0(:, i), n(:, i));
            t0 = LineLineIntersect(v0(:, j), n0(:, j), v0(:, i), n(:, i));
            t1 = LineLineIntersect(v1(:, j), n1(:, j), v0(:, i), n(:, i));
            
            % build occluding interval
            if (0.0 < t0(1)) && (0.0 < t1(1))
                % ++ no tweaks
            elseif (-d0(j) <= t0(1)) && (t0(1) <= 0.0) ...
                && (-d1(j) <= t1(1)) && (t1(1) <= 0.0)
                % 00
                continue
            elseif (t0(1) < -d0(j)) && (t1(1) < -d1(j))
                % --
                continue
            elseif (0.0 < t0(1)) && (-d1(j) <= t1(1)) && (t1(1) <= 0.0)
                % +0
                t1(2) = t(2);
            elseif (-d0(j) <= t0(1)) && (t0(1) <= 0.0) && (0.0 < t1(1))
                % 0+
                t0(2) = t(2);
            elseif (0.0 < t0(1)) && (t1(1) < -d1(j))
                % +-
                if dot(nn(:, j), n(:, i)) >= 0.0
                    t1(2) = inf();
                else
                    t1(2) = -inf();
                end
            elseif (t0(1) < -d0(j)) && (0.0 < t1(1))
                % -+
                if dot(nn(:, j), n(:, i)) >= 0.0
                    t0(2) = inf();
                else
                    t0(2) = -inf();
                end
            elseif (-d0(j) <= t0(1)) && (t0(1) <= 0.0) && (t1(1) < -d1(j))
                % 0-
                t0(2) = t(2);
                if dot(nn(:, j), n(:, i)) >= 0.0
                    t1(2) = inf();
                else
                    t1(2) = -inf();
                end
            elseif (t0(1) < -d0(j)) && (-d1(j) <= t1(1)) && (t1(1) <= 0.0)
                % -0
                t0(2) = t(2);
                if dot(nn(:, j), n(:, i)) >= 0.0
                    t1(2) = inf();
                else
                    t1(2) = -inf();
                end
            else
                warning('Impossible!?');
            end
            
            % modify visible interval
            interval = IntervalDifference(interval, t0(2), t1(2));
            
        end
        
        % set inteval
        visibility{i} = interval;
        
%         % plot interval
%         for j = 1 : size(interval, 2)
%             x0 = n(:, i) * interval(1, j) + v0(:, i);
%             x1 = n(:, i) * interval(2, j) + v0(:, i);
%             hold on
%             line([x0(1), x1(1)], [x0(2), x1(2)], 'LineWidth', 8, 'Color', 'r');
%         end
        
    end
    
end




function TxVisibility(point, tx, geometry)
%     % an alternative line formulation, may offer more elegant solution
%     % uses a normal to the line. also f(x, y) gives signed distance to (x,y)
%     % f(x, y) = A * x + B * y + C
%     v0 = verts(:, edges(1, :));
%     v1 = verts(:, edges(2, :));
%     A = v0(2, :) - v1(2, :);
%     B = v1(1, :) - v0(1, :);
%     C = v0(1, :) .* v1(2, :) - v1(1, :) .* v0(2, :);
%     D = (A .^ 2 + B .^ 2) ^ 0.5;
%     A = A ./ D;
%     B = B ./ D;
%     C = C ./ D;
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


function [] = DebugPlot(tx, rx, geometry, depth, rays)
    
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
    nrays = numel(rays);
    for i = 1 : nrays
        line(ah, 'XData', rays{i}(1, :), 'YData', rays{i}(2, :), 'Color', 'r');
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
