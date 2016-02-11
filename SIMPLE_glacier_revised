% 1-D staggered-grid glacier modelwith cole
% written june 12 2015 rsa and cole
% all units in SI units

clearvars
figure(1)
clf
figure(2)
clf

%% initialize
rho_i = 917;
g     = 9.81;
A     = 2.1e-16; % Pa-3 yr-1

% initial topography
zbmax = 4000;
zbmin = 2000;
xstar = 7000;

dx = 50;
xmax = 20000;
x = dx/2:dx:xmax-(dx/2);
xedge = 0:dx:xmax;

zb = zbmin+(zbmax-zbmin)*exp(-x/xstar);
H = zeros(size(x));
%Hedge = zeros(size(xedge));
z = zb+H;

% now set up meteorology
ELA = 2900;
dbdz = 0.01;
bcap = 2.0; % m/yr
b = dbdz*(z-ELA);
b = min(b,bcap);

dt = 0.001;
tmax = 400;

t = 0:dt:tmax;
imax = length(t);
nplots = 20;
tplot = tmax/nplots;

nframe = 0;

%% run

for i = 1:imax
    
b = dbdz*(z-ELA); % local net balance calculated at cell centers at ice surface
b = min(b,bcap);
Hedge = H(1:end-1)+0.5*diff(H); % interpolates ice thickness to cell edges
S = abs(diff(z)/dx); % slope of ice surface calculated at cell edges

Q = (A/5).*((rho_i*g*S).^3).*(Hedge.^5); % ice discharge due to internal deformation
Q = [0 Q 0]; %takes care of boundary conditions

dHdt = b - (diff(Q)/dx); % ice continuity
H = H + (dHdt*dt); %updates ice thickness
H = max(H,0);

z = zb+H; %updates topography

% now for some plotting
if rem(t(i),tplot)==0
    nframe=nframe+1;

    figure(1)
    %subplot(1,2,1)
    subplot('position',[0.1 0.1 0.7 0.85])
        %subplot('position',[left bottom width height])
    plot(x/1000,zb,'k','linewidth',2)
    hold on
    plot(x/1000,z)
    plot(x/1000,ELA*ones(size(x)),'g--','linewidth',2)
    axis([0 xmax/1000 zbmin zbmax+200])
    xlabel('Horizontal Distance (km)','fontname','arial','fontsize',18)
    ylabel('Elevation (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    
    %subplot(1,2,2)
    subplot('position',[0.85 0.1 0.1 0.85])
    plot(b,z,'b','linewidth',2)
    hold on
    plot(zeros(size(z)),z,'g--')
    axis([min(b) 1.5*bcap zbmin zbmax+200])
    xlabel('b (m/yr)','fontname','arial','fontsize',18)
    %ylabel('Elevation (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    
    Qanal = cumsum(b)*dx; %analytic solution for ss ice discharge
    Qanal = max(Qanal,0);
    
    figure(2)
    plot(xedge/1000,Q/1000)
    hold on
    plot(x/1000,Qanal/1000,'g--')
    axis([0 xmax/1000 0 bcap*xmax/2000])
    xlabel('Horizontal Distance (km)','fontname','arial','fontsize',18)
    ylabel('Ice discharge (1000 m^2/yr)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    pause(0.2)
        
end

end



%% finalize



