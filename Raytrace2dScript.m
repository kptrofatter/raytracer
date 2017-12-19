% %============================================================================%
% % Duke University                                                            %
% % K. P. Trofatter                                                            %
% % kpt2@duke.edu                                                              %
% %============================================================================%
% Raytrace2dScript() - runs a raytrace simulation.

%% Environment
close('all');
clear();
clc();

% physical constants
ep0 = 8.85e-12;        % [F/m]
mu0 = pi * 4.0e-7;     % [H/m]
c0  = 1.0 / sqrt(ep0 * mu0); % [m/s]


%% Parameters
f0 = 20.0e9;              % [Hz] frequency of operation
w0 = 2.0 * pi * f0;       % [rad/s] angular frequency
beta = w0 / c0;           % [rad/m] wavenumber
lambda = 2.0 * pi / beta; % [m] wavelength
alpha = 1.1;              % [Np/m] attenuation
E0 = 1.0;                 % [V/m] input field

% antenna gain patterns
rng(0);
G_tx = 1.0 * exp(1.0j * 2.0 * pi * zeros(1, 32)); % [#] tx complex gain
G_rx = 1.0 * exp(1.0j * 2.0 * pi * zeros(1, 32)); % [#] tx complex gain

% experiment
switch 'waveguide'
case 'dreza'
    rx = [1.0; 0.1]; % [m] rx position
    tx = [5.4; 3.0]; % [m] tx position
    verts_scene = [ ... % [m] scene vertices
        0.0, 3.5, 3.5, 5.5, 5.5, 0.5, 0.5, 0.0; ...
        0.0, 0.0, 2.0, 0.0, 4.0, 4.0, 3.5, 3.5];
    edges_scene = [ ... % [m] scene edge vertex indices
        1, 2, 2, 4, 5, 6, 7, 8; ...
        2, 3, 4, 5, 6, 7, 8, 1];
    npoints = 8;
    theta = linspace(0.0, 2.0 * pi, npoints + 1);
    theta(end) = [];
    verts_object = [ ... % [m] object vertices
        cos(theta); ...
        sin(theta)];
    edges_object = [ ... % [m] object edges
        1, 2 : npoints; ...
        2 : npoints, 1];
    depth = 4; % [#] max recursion depth
    
case 'waveguide'
    rx = [0.0; 0.09]; % [m] rx position
    tx = [1.0; 0.0]; % [m] tx position
    verts_scene = [ ... % [m] scene vertices
        0.0,  0.0,  1.0, 1.0; ...
        0.1, -0.1, -0.1, 0.1];
    edges_scene = [ ... % [m] scene edges
        1, 2; ...
        4, 3];
    verts_object = [-1.0, 1.0; 0.0, 0.0]; % [m] object vertices
    edges_object = [1; 2]; % [m] object edges
    depth = 2; % [#] max recursion depth

case 'bend'
    rx = [0.0; 0.0]; % [m] rx position
    tx = [1.4; -0.4]; % [m] tx position (PUT y=[2.0,-0.2] for error)
    verts_scene = [ ... % [m] scene vertices
        0.0,  0.0,  0.9, 1.0, 1.4, 1.6; ...
        0.1, -0.1, -0.1, 0.1, -0.6, -0.4];
    edges_scene = [ ... % [m] scene edges
        1, 2, 4, 3; ...
        4, 3, 6, 5];
    verts_object = zeros(2, 0); % [m] object vertices
    edges_object = zeros(2, 0); % [m] object edges
    depth = 32; % [#] max recursion depth

case 'retro'
    rx = [2.0; 6.0]; % [m] rx position
    tx = [3.0; 5.0]; % [m] tx position
    verts_scene = [ ... % [m] scene vertices
        0.0, 8.0, 0.0; ...
        0.0, 0.0, 8.0];
    edges_scene = [ ... % [m] scene edges
        1, 1; ...
        2, 3];
    verts_object = zeros(2, 0); % [m] object vertices
    edges_object = zeros(2, 0); % [m] object edges
    depth = 2; % [#] max recursion depth

end
    verts_object = zeros(2, 0); % [m] object vertices
    edges_object = zeros(2, 0); % [m] object edges
% debug ray trace plot
debug = 'distance'; % debug plot flag

% timeline
nframes = 100; % [#]
t0 = 0.0;    % [s] start time
t1 = 1.0;    % [s] end time
dt = t1 - t0;
tstep = dt / nframes; % [s] time step
T  = linspace(t0, t1, nframes); % time vector


%% Run
% initiate rx field
E = zeros(1, nframes);

% simulation loop
for i = 1 : numel(T)
    
    % progress
    fprintf('[%i/%i] Frame\n', i, numel(T));
    
    % animate object
    S = 0.05 * [1.0, 0.0; 0.0, 1.0];
    theta = 2.0 * pi * (T(i) - t0) / dt / 2;
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    A.M = R * S;
    %A.v = [4.5; 1.5]; % dreza
    A.v = [0.5; 0.0]; % waveguide
    obj = AffineXform(A, verts_object);
    
    % assemble geometry
    offset = size(verts_scene, 2);
    geometry.verts = [verts_scene, obj];
    geometry.edges = [edges_scene, edges_object + offset];
    
    % ray trace
    traces = Raytrace2d(rx, tx, geometry, depth, debug);
    drawnow();
    
    % compute trace angles at tx and rx
    ntraces = numel(traces);
    theta_tx = zeros(1, ntraces);
    theta_rx = zeros(1, ntraces);
    for j = 1 : ntraces
        
        % get trace
        trace = traces{j};
        
        % tx exit angle
        dx = trace(:, end - 1) - trace(:, end);
        theta_tx(j) = atan2(dx(2), dx(1));
        % rx entrance angle
        dx = trace(:, 1) - trace(:, 2);
        theta_rx(j) = atan2(dx(2), dx(1));
        
    end
    
    % construct histogram of tx angles
    edges_tx = linspace(-pi, pi, 2 * numel(G_tx) + 1);
    edges_tx = edges_tx(1 : 2 : end);
    bins_tx = histcounts(theta_tx, edges_tx);
    % construct histogram of rx angles
    edges_rx = linspace(-pi, pi, 2 * numel(G_rx) + 1);
    edges_rx = edges_rx(1 : 2 : end);
    bins_rx = histcounts(theta_rx, edges_rx);
    
    % interfere wave
    for j = 1 : ntraces
        
        % get trace
        trace = traces{j};
        
        % compute distance
        d = 0.0;
        for k = 1 : size(trace, 2) - 1
            d = d + norm(trace(:, k) - trace(:, k + 1));
        end
        
        % bin angle at tx
        dx = trace(:, end - 1) - trace(:, end);
        a = atan2(dx(2), dx(1));
        ibin_tx = (edges_tx(1 : end - 1) <= a) & (a <= edges_tx(2 : end));
        % bin angle at rx
        dx = trace(:, 1) - trace(:, 2);
        a = atan2(dx(2), dx(1));
        ibin_rx = (edges_rx(1 : end - 1) <= a) & (a <= edges_rx(2 : end));
        
        % propagate trace
        k = beta - 1.0j * alpha;
        Et = E0 * (G_tx(ibin_tx) / bins_tx(ibin_tx));
        Et = Et * exp(-1.0j * k * d);
        Et = Et / (G_rx(ibin_rx) / bins_rx(ibin_rx));
        
        % accumulate
        E(i) = E(i) + Et;
        
    end
    
    % animate
    %Gif(fh, gif, delay, loopcount);
end


%% Draw
% draw E(t)
fh = figure(2);
clf(fh);
% abs(E)
ah = subplot(1, 2, 1, 'Parent', fh);
plot(ah, T, abs(E));
title(ah, '|E(t)| @ rx');
xlabel(ah, 't[s]');
ylabel(ah, '|E(t)|[V/m]');
grid(ah, 'on');
xlim(ah, [T(1), T(end)]);
% angle(E)
ah = subplot(1, 2, 2, 'Parent', fh);
plot(ah, T, angle(E));
title(ah, 'arg(E(t)) @ rx');
xlabel(ah, 't[s]');
ylabel(ah, 'arg(E(t))[rad]');
grid(ah, 'on');
xlim(ah, [T(1), T(end)]);

% draw tx gain pattern
fh = figure(3);
clf(fh);
% draw modulus of G_tx
ah = subplot(1, 2, 1, 'Parent', fh);
theta = linspace(0.0, 2.0 * pi, numel(G_tx) + 1);
theta(end) = [];
plot(theta, abs(G_tx));
title(ah, '|G_{tx}(\theta)|');
xlabel(ah, '\theta[rad]'); ylabel(ah, '|G_{tx}|');
grid(ah, 'on');
xlim(ah, [T(1), T(end)]);
% draw argument of G_tx
ah = subplot(1, 2, 2, 'Parent', fh);
plot(theta, angle(G_tx));
title(ah, 'arg(G_{tx}(\theta))');
xlabel(ah, '\theta[rad]'); ylabel(ah, 'arg(G_{tx})[rad]');
grid(ah, 'on');
xlim(ah, [T(1), T(end)]);

% draw rx gain pattern
fh = figure(4);
clf(fh);
% draw modulus of G_rx
ah = subplot(1, 2, 1, 'Parent', fh);
theta = linspace(0.0, 2.0 * pi, numel(G_rx) + 1);
theta(end) = [];
plot(theta, abs(G_rx));
title(ah, '|G_{rx}(\theta)|');
xlabel(ah, '\theta[rad]'); ylabel(ah, '|G_{rx}|');
grid(ah, 'on');
xlim(ah, [T(1), T(end)]);
% draw argument of G_rx
ah = subplot(1, 2, 2, 'Parent', fh);
plot(theta, angle(G_rx));
title(ah, 'arg(G_{rx}(\theta))');
xlabel(ah, '\theta[rad]'); ylabel(ah, 'arg(G_{rx})[rad]');
grid(ah, 'on');
xlim(ah, [T(1), T(end)]);

% 2016a or later
% % draw polar of modular radiation pattern
% theta = linspace(0.0, 2.0 * pi, numel(G_tx) + 1);
% theta(end) = [];
% polarplot(theta, abs(G_tx));
% title(ah, 'arg(G_tx)');
% 
% % draw polar of modular radiation pattern
% theta = linspace(0.0, 2.0 * pi, numel(G_rx) + 1);
% theta(end) = [];
% polarplot(theta, abs(G_rx));
% title(ah, 'arg(G_rx)');


%==============================================================================%
%                                                                              %
%                                                                              %
%                                                                              %
%==============================================================================%
