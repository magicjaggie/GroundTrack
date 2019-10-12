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