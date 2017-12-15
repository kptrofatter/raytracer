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
ep0 = 8.85e-12;           % [F/m]
mu0 = pi*4.0e-7;          % [H/m]
c0  = (ep0 * mu0) ^ -0.5; % [m/s]

%% Parameters
f0 = 20.0e9;             % [Hz] frequency of operation
w0 = 2.0 * pi * f0;      % [rad/s] angular frequency
k0 = w0 / c0;            % [rad/m] wavenumber
lambda0 = 2.0 * pi / k0; % [m] wavelength

%rx = [1.0; 0.0]; % [m] rx position
rx =  [1.0; 1.5];

tx = [5.5; 3.0]; % [m] tx position
tx = [1.2; 0.5];

geometry.verts = [ ... % [m] scene vertices
    0.0, 3.5, 3.5, 5.5, 5.5, 0.5, 0.5, 0.0,   1.0 2.0 1.0 2.0; ...
    0.0, 0.0, 2.0, 0.0, 4.0, 4.0, 3.5, 3.5,   1.0 1.0 2.0 2.0];

geometry.edges = [ ... % [m] scene edges
    1, 2, 2, 4, 5, 6, 7, 8,  9 ,10 ; ...
    2, 3, 4, 5, 6, 7, 8, 1,  12,11  ];

depth = 4; % [#] max recursion depth

debug = true(); % debug plot flag

%% Run

traces = Raytrace2d(rx, tx, geometry, depth, debug);


% ray trace
for i = -1.6 : 0.04 : 3.6
    traces = Raytrace2d([i + 2; i], tx, geometry, depth, debug);
    drawnow();
end

% interfere rays
    % tx gain
    % tx phase
    % phase propagation
    % magnitude attenutation
    % rx gain
    
% plot pattern

%==============================================================================%
%                                                                              %
%                                                                              %
%                                                                              %
%==============================================================================%
