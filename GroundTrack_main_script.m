clc;
clear;

% Nb of revolutions
N = 3.65;

% Input parameters
a = 8350;
e = 0.19760;
i = 60;
Omega = 270; % deg
omega = 45; % deg
f_0 = 230; % deg - theta
omega_E = 15.04; % deg/h - Rotation of the Earth
mu = 398600.44; % orbital period of the satellite

% deg/s
omega_E = omega_E / 3600;

% Getting the cartesian coordinates
[r,v] = kep2car(a, e, i, Omega, omega, f_0, mu);

state = [r;v];



%% Calculating T and Y using Lab1

% Time
Tperiod = 2*pi()*sqrt(a^3/mu); % period, s
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
lambda0 = -180;

[alpha, delta, lon, lat] = groundTrack(Y, lambda0, Tperiod, omega_E, N);
lon = flip(lon);
lat = flip(lat);

% Plot the ground tracks

figure(1)

image_file = '/Users/morgane/Desktop/Laboratories Orb Mech/GroundTrack/Earth.jpg';
% img = imread(image_file);
% image('CData', img, 'XData', [-180 180], 'YData', [-90 90]);
% %imshow(img)

% geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
% title('Satellite Ground Track')

cdata = (imread(image_file));
imagesc([-180,180],[-90, 90],cdata);

hold on;
L = length(Y(:,1));

% Plot the complete revolutions
for k = 0:floor(N)-1
    plot(lon(1 +k*L : L +k*L),lat(1 +k*L : L +k*L),'-r');
end

% Plot the last partial revolution
if N>floor(N)
    d = N - floor(N); % Percentage of the last revolution that is done
    plot(lon(1 +floor(N)*L : floor(L*d) +floor(N)*L),lat(1 +floor(N)*L : floor(L*d) +floor(N)*L), '-r');
end
    
xlabel("Longitude");
ylabel("Latitude");