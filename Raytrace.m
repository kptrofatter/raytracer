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
    
function [rays] = Raytrace(tx, rx, geometry, depth, debug)
    
    rays = [];
    
    if debug
        DebugPlot(tx, rx, geometry, depth, rays)
    end
    
end


function [t] = EdgeEdgeIntersect(c1, n1, c2, n2)
    
    % parallel test
    if (abs(dot(n1, n2)) - 1.0) < 1.0e-8
        % parallel, no intersect, not even if overlapping
        t = inf(2, 1);
    else
        % compute intersection along each line
        t = [n1, n2] \ (c2 - c1);
    end
    
end


function [visibility] = Visibility(point, geometry)
    
    % get geometry
    verts = geometry.verts;
    edges = geometry.edges;
    
    % count geometry
    nverts = size(verts, 2);
    nedges = size(edges, 2);
    
    % compute edge line equations
    % x(t) = t * n + v0
    v0 = verts(:, edges(1, :)); % edge origin
    v1 = verts(:, edges(2, :));
    n = v1 - v0;
    d = sum(n .^ 2) .^ 0.5; % edge length
    n = n ./ repmat(d, [2, 1]); % edge direction
    
    % compute point-to-edge-ends line equations
    P = repmat(point, [1, nverts]);
    % end 0
    n0 = v0 - P;
    d0 = sum(n0 .^ 2) .^ 0.5;
    n0 = n0 ./ repmat(d0, [2, 1]);
    % end 1
    n1 = v1 - P;
    d1 = sum(n1 .^ 2) .^ 0.5;
    n1 = n1 ./ repmat(d1, [2, 1]);
    
    % initiate visibility
    visibility = cell(1, nedges);
    for i = 1 : nedges
        visibility{i} = Interval(0.0, dx(i));
    end
    
    % edge visibilty loop
    for i = 1 : nedges
        
        % get edge interval
        interval = visibility{i};
        
        % occluding edges loop
        for j = 1 : nedges
            
            % skip self-test
            if i == j
                continue
            end
            
            % test for empty interval
            if isempty(interval)
                break
            end
            
            % occlude
            %interval = Occlude(interval, point, i, j, geometry);
            

            EdgeEdgeIntersect(c1, n1, c2, n2)
        end
        
        % set inteval
        visibility{i} = interval;
        
    end
    
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
        elseif (a < c) && (d < b)
            % a cd b
            intersect(:, idiff) = [c; d];
            idiff = idiff + 1;
        elseif (c < a) && (b < d)
            % c ab d
            intersect(:, idiff) = [a; b];
            idiff = idiff + 1;
        elseif (a < c) && (b < d)
            % a c b d
            intersect(:, idiff) = [c; b];
            idiff = idiff + 1;
        elseif (c < a) && (d < b)
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
        if (b < c) || (d < a)
            % ab cd, cd ab
            diff(:, idiff) = [c; d];
            idiff = idiff + 1;
        elseif (a < c) && (d < b)
            % a cd b
        elseif (c < a) && (b < d)
            % c ab d
            diff(:, idiff + 0) = [c; a];
            diff(:, idiff + 1) = [b; d];
            idiff = idiff + 2;
        elseif (a < c) && (b < d)
            % a c b d
            diff(:, idiff) = [b; d];
            idiff = idiff + 1;
        elseif (c < a) && (d < b)
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
    fh = figure();
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
    line(ah, 'Xdata', e(1, :), 'Ydata', e(2, :), ...
        'Color', 'k', 'LineWidth', 2.0);
    
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
