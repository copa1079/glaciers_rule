% 1d glacier modelwith cole
% written june 12 2015 rsa and cole
% alll units in SI units

clear all
figure(1)
clf
figure(2)
clf

%% initialize
rho_i = 917; %kgm^3
g = 9.81; %ms^-2
A = 2.1e-16; % Pa^-3 yr^-1

% initial topography (meters)
zbmax = 4000;
zbmin = 2000;
xstar = 10000;

dx = 100;
xmax = 20000;
x = dx/2:dx:xmax-(dx/2);
xedge = 0:dx:xmax;

zb = zbmin+(zbmax-zbmin)*exp(-x/xstar);
H = zeros(size(x));
%Hedge = zeros(size(xedge));
z = zb+H;

figure(1)
plot(x/1000,zb)

% now set up meteorology
ELA = 3200; % meters
dbdz = 0.01; %T^-1
bcap = 2; % m/yr
b = dbdz*(z-ELA); % m/yr
b = min(bcap,b); % m/yr

dt = 0.001;
tmax = 500;

t = 0:dt:tmax;
imax = length(t);
nplots = 20;
tplot = tmax/nplots;

nframe = 0;

%% run

for i = 1:imax
    
b = dbdz*(z-ELA);
Hedge = H(1:end-1)+0.5*diff(H);

S = abs(diff(z)/dx); % slope of ice surface
Q = (A/5).*((rho_i*g*S).^3).*(Hedge.^5);

Q =[0 Q 0]; %takes care of boundary conditions
dHdt = b - (diff(Q)/dx);

H = H + (dHdt*dt);
H = max(H,0);

z = zb+H;

if rem(t(i),tplot)==0
    nframe = nframe + 1
    figure(1)
    plot(x/1000,zb,'k','linewidth',2)
    hold on
    plot(x/1000,z,'b')
    plot(x/1000,ELA*ones(size(x)),'g--')
    axis([0 xmax/1000 zbmin zbmax+200])
    xlabel('Horizontal Distance (km)','fontname','arial','fontsize',18)
    ylabel('Elevation (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    
    figure(2)
    plot(xedge/1000,Q) %discharge curve set
    hold on
    plot(x/1000,cumsum(b)*dx,'g--') %mass balance curve set1
    axis([0 18 0 40000])
    xlabel('Horizontal Distance (km)','fontname','arial','fontsize',18)
    ylabel('Discharge (km^2/s)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    pause(0.5)
    
    %figure(3)
    %plot(z/1000,b,'r')
    %xlabel('Mass balance (m)','fontname','arial','fontsize',18)
    %ylabel('Vertical Distance (km)','fontname','arial','fontsize',18)
    %set(gca,'fontsize',14,'fontname','arial')

end

end



%% finalize

figure(2)
plot(xedge/1000,Q)


