clc;
clear;

% Nb of revolutions
k = 3.25;

% Input parameters
a = 8350;
e = 0.19760;
i = 60;
Omega = 270; % deg
omega = 45; % deg
f_0 = 230; % deg - theta
omega_E = 15.04; % deg/h - Rotation of the Earth
mu = 398600.44; % orbital period of the satellite


% Getting the cartesian coordinates
[r,v] = kep2cart(a, e, i, Omega, omega, f_0, mu);

state = [r;v];



%% Calculating T and Y using Lab1

% Time
Tperiod = 2*pi()*sqrt(a^3/mu); % period
tfin = Tperiod;
t0 = 0;
tspan = linspace(0,Tperiod, 100)'; % faster
%tspan = [t0,k*tfin]; % slower

% Coordinates
Y0 = [r(1) r(2) r(3) v(1) v(2) v(3)];

% Set options
options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-12 );

% Perform the integration
[ T, Y ] = ode113( @(t,y)ode_keplerian_orbit(t,y,mu), tspan, Y0, options);
r = [Y(:,1),Y(:,2),Y(:,3)];
v = [Y(:,4),Y(:,5),Y(:,6)];


% Initial condition
lambda0 = 0;

[alpha, delta, lon, lat] = groundTrack(Y, lambda0, Tperiod, omega_E);
