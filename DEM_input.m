%% SCRIPT TO GENERATE BEDROCK TOPOGRAPHY FOR USE IN 2D ARKANSAS GLACIER MODEL

%  Cole C. Pazar, November 16th 2015

clear all
close all

%% SET TOGGLE TO CHOOSE INPUT MODE
toggle = 1; % 1 = synthetic topography, 2 = real DEM

if toggle == 1
    %% PARABOLIC (Y) EXPONENTIAL (X) VALLEY TEST CASE
    % Generating x-y plane
    xMax = 1e4;
    dx = 100;
    yMax = 2.5e3;
    dy = 100;

    y = -yMax:dy:yMax;
    x = 0:dx:xMax;

    [X Y] = meshgrid(x,y);
    Z = zeros(length(y),length(x));

    % Generating z elevations
    zMax = 1e3; % divide elevation
    xStar = 5e3; % characteristic lengthscale along x
    clRow = find(y==0); % index for glacier centerline
    power = 2; % power for cross-valley profile
    damper = 2e-4; % squish-factor for cross-valley profile
    Z = damper*Y.^power;
    Z = Z + zMax*exp(-X/xStar); % making centerline profile


    figure(1)
    clf
    %orient landscape
    subplot(2,2,1:2)
    surf(X,Y,Z)
    axis equal
    xlabel('Down-valley distance [m]','fontsize',16)
    ylabel('Cross-valley distance [m]','fontsize',16)
    zlabel('Vertical distance [m]','fontsize',16)
    set(gca,'fontsize',14)

    subplot(2,2,3)
    plot(x,Z(clRow,:),'linewidth',3)
    xlabel('Down-valley distance [m]','fontsize',18)
    ylabel('Vertical distance [m]','fontsize',18)
    set(gca,'fontsize',14)
    axis equal

    subplot(2,2,4)
    plot(Y(:,1),Z(:,1),'b','linewidth',3)
    hold on
    plot(Y(:,1),Z(:,ceil(length(Z)/2)),'g','linewidth',3)
    plot(Y(:,1),Z(:,length(Z)),'r','linewidth',3)
    xlabel('Cross-valley distance [m]','fontsize',18)
    ylabel('Vertical distance [m]','fontsize',18)
    legend('x = 0','x = xMax/2','x = xMax','location','north')
    set(gca,'fontsize',14)
    axis equal

%% READ IN THE REAL DEM
elseif toggle == 2
    z = imread('/Users/ccp127/Desktop/Thesis/Figures/DEM.tiff'); 
%%%%%%%%%%% this is the part of the code where you choose the location of
%%%%%%%%%%% the DEM file in .tiff format. 
    pxSize = 26.33239;
    numEpx = size(z,1);
    numNpx = size(z,2);

    figure(1)
    orient landscape
    h = imagesc(0:pxSize:numNpx*pxSize,0:pxSize:numEpx*pxSize,z);
    f = colorbar;
    set(f,'ylim',[2900 4300])
    set(gca,'YDir','reverse','fontsize',14)
    set(h,'AlphaData',z~=0)
    axis image
    ylabel(f,'Elevation [m]','fontsize',18)
    xlabel('Relative easting [m]','fontsize',18)
    ylabel('Relative northing [m]','fontsize',18)
end
