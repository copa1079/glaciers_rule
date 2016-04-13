clear global
clearvars
% author: Cole C. Pazar, April 3rd, 2016
%% BEGIN CODE: ORGANIZE ALL CONSTANTS AND VALUES [=] SI units
S = 0.01; % approx. topographic slope [=] m/m derived from QGIS profile tool
g = 9.81; % gravity
Cd = 0.4; % approx. drag coefficent for a rough boulder in high-Re
k  = 0.4;  % von karman's constant
rhorock = 2650;  % density of granite
rhowater = 1000; % denisty of water
    deltarho = rhorock-rhowater; % bouyant density of the boulders
width1 = 80;     % channel width approximation [=] meters
width2 = width1; % channel width approximation [=] meters
D84 = 0.75;      % 84th percentile of the grain sizes on the bed
Ro = D84/2;      % half that for the radius
z_o = 0.1*D84;   % bed roughness

%% BOULDER DATA
N1 = 47; % number of boulders measured on terrace 1
N2 = 24; % number of boulders measured on terrace 1
D1 = 3.75; % mean boulder diameter for terrace 1 (statsitical average)
D2 = 6.80; % mean boulder diameter for terrace 2 (statsitical average)
    R1 = D1/2; % mean boulder radius for terrace 1
    R2 = D2/2; % mean boulder radius for terrace 2
    rad1 = asin(Ro/(R1+Ro)); % angle1 in geometry of torque [=] radians
    rad2 = asin(Ro/(R2+Ro)); % angle2 in geometry of torque [=] radians
    angle1 = rad1*(180/pi); 
    angle2 = rad2*(180/pi);
    cosangle1 = cos(rad1);
    cosangle2 = cos(rad2);
    
%% UNCERTAINTIES IN FLOOD BOULDER DATA
    uD1 = 0.24; % uncertainty in D1 --> standard deviation of the mean
    uD2 = 0.33; % uncertainty in D2 --> standard deviation of the mean
    uR1 = uD1/2; % uncertainty in R1
    uR2 = uD2/2; % uncertainty in R2
    r1 = (R1*Ro)/(R1+Ro); % gravtitationl torque lever arm 1
    r2 = (R2*Ro)/(R2+Ro); % gravtitationl torque lever arm 2

%% FLOODWATERS HEIGHT CALCULATIONS:
height1 = (1/S)*(k^2)*(1/(Cd*1.2))*(Ro/(R1+Ro))*(D1/cosangle1)*(deltarho/rhowater);
height2 = (1/S)*(k^2)*(1/(Cd*1.2))*(Ro/(R2+Ro))*(D2/cosangle2)*(deltarho/rhowater);

%% PALEOLAKE GEOMETRIES DERIVED FROM FIELDWORK, GOOGLE EARTH, AND QGIS:

A1 = 7.4; % area of clear creek lake [=] km^2
    h1 = 0.025; % mean depth of paleo clear creek lake [=] km
A2 = 45; % area of Clear Creek lake [=] km^2
    h2 = 0.04;  % mean depth of the Granite lake [=] km
v1 = A1*h1; % approximate volume of the paleo Clear Creek lake [=] km^3
    V1 = v1*10^9; % converting to m^3
v2 = A2*h2; % approximate volume of the Granite lake [=] km^3
    V2 = v2*10^9; % converting to m^3

%% FINALIZE
velocity1 = (1/0.4) * (sqrt(9.81*height1*S)) * (log(height1/z_o)-1);
velocity2 = (1/0.4) * (sqrt(9.81*height1*S)) * (log(height1/z_o)-1);
Q1 = width1 * height1 * velocity1; % discharge at ~17 ka
Q2 = width2 * height2 * velocity2; % discharge at ~19 ka
drain_time1 = (V2/Q1)/3600; % in hours
drain_time2 = (V2/Q2)/3600; % in hours

%% PRINT LINE
%                                   output:
% depths       
height1 %                           = 27.8901 meters  
height2 %                           = 29.8696 meters  
% discharges
Q1 %                                = 4.5381e+04 cubic meters per second 
Q2 %                                = 4.8602e+04 cubic meters per second
% time to drain the 'Granite' lake
drain_time1 %                       = 11.0178 hours       
drain_time2 %                       = 10.2876 hours
%
% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
