%% DATE: Feb. 10th 2016
%% CLEANED UP GLACIER MODEL FOR CLEAR CREEK
%  1D FTCS STAGGERED GRID NUMERICAL MODEL: CLEAR CREEK GLACIER
%  COLE C. PAZAR AND ROBERT S. ANDERSON
clearvars
figure(1)
clf
figure(2)
clf
% figure(3)
% clf

%% initialize
step = 200;
rho_i = 917;
g = 9.81;
A = 2.1e-16; % Pa-3 yr-1
slide_ratio = 0.5;% ratio of sliding speed to internal defm speed

xmax = 30000; % was 28 km ...
dx = xmax/step;
x = dx/2:dx:xmax-(dx/2);
xedge = 0:dx:xmax;

% now the Clear Creek, Colorado case
zb0 = 3425;
slope0 = 0.0285; % was 0.0338
zb1 = 516;
xstar = 900; % was 1687
zbline = zb0 - (slope0*x);
zbexp = zb1 * exp(-x/xstar);

zb = zbline + zbexp;

% valley width as a function of distance downvalley
% turn this off by setting Wmin=W0;
W0 = 3000;
Wstar = 4000;
Wmin = 1700;
W = Wmin + (W0-Wmin)*exp(-x./Wstar);
Wedge = W(1:end-1)+0.5*diff(W); % interpolates valley width to cell edges
Wedge = [Wedge(1) Wedge Wedge(end)];

% now put in a wall at the bottom of the valley akin to the arkansas case
% in clear ck and pine creek glaciers
slopewall = 0.3;
xwall = 29450;
zbwall = 2800 + (slopewall * (x-xwall));
zb = max(zb,zbwall);
beyond = find(x>29450);
slope2 = 0.01;
zb(beyond) = 2900+(slope2*(x(beyond)-xwall));

dambase = find(zb==min(zb));
xdambase = x(dambase);
zbmin = min(zb);
zbmax = max(zb);

H = zeros(size(x));
z = zb+H;

% now set up meteorology (some based upon brugger's suggestion in sawatch
% paper 2006)
ELA0 = 3300;      % AVERAGE ELA: 3270,   SEE BRUGGER (2010)
sigma_ELA = 300;  % used to be 250

dbdz = 0.006; % m/y/m, typically 0.01
bcap = 1.25;  % m/yr, usually 1.25-2.00
b0 = dbdz*(z-ELA0);
b0 = min(b0,bcap);
minzb = find(zb==min(zb));
minb = b0(minzb);

% time array
dt = 0.0025;
tmax = 4000; % max time interval of growth, 
t = 0:dt:tmax;

% plotting controls 
imax = length(t);
nplots = 80;
tplot = tmax/nplots;
nframe = 0;

% control the Climate:
big_period = 100000; % years for the 100 ka milankovich cycle
med_period = 2500; % years so that the peak is centered
sol_period = 11; % years for solar cycle

big_shift = -50000; % years
med_shift = 1000;

sol_shift = 11*randn(size(t));  % include a solar cycle  w/  1 m / yr

randomsize_t = randn(size(t));

big = (sigma_ELA/4)*sin(2*pi*(t+big_shift)/big_period);
medium = (sigma_ELA/4).*sin(2*pi*(t+med_shift)/med_period);
solar = (sol_shift).*sin(2*pi*t/sol_period);

ELA = ELA0 + sigma_ELA + big + medium + solar - 180;

%% run

for i = 1:imax
    
b = dbdz*(z-ELA(i)); % local net balance calculated at cell centers at ice surface
b = min(b,bcap);

Hedge = H(1:end-1)+0.5*diff(H); % interpolates ice thickness to cell edges
S = abs(diff(z)/dx); % slope of ice surface calculated at cell edges

Udef = (A/5).*((rho_i*g*S).^3).*(Hedge.^4); %mean defm speed
Q = (A/5).*((rho_i*g*S).^3).*(Hedge.^5); % ice discharge due to internal deformation
Qsl = slide_ratio * Udef.*Hedge;
Q = Q + Qsl;
Q =[0 Q 0]; %takes care of boundary conditions

%dHdt = b - (diff(Q)/dx); % ice continuity
dHdt = b - (1./W).*(diff(Q.*Wedge)/dx); %continuity allowing width to vary
H = H + (dHdt*dt); %updates ice thickness
H = max(H,0);

z = zb+H; %updates topography
glacier = find(H>0);
term(i)=x(glacier(end));

% track dam height
if(x(glacier(end))>xdambase)
    dam_ht(i)=zb(glacier(end))-min(zb);
else
    dam_ht(i)=0;
end
    dam_ht(i)=max(0,dam_ht(i));
    
% now for some plotting
if rem(t(i),tplot)==0
    nframe=nframe+1;
    
    M(:,nframe) = getframe(gcf); 
    
    load PK_profile.txt
    zb_ccp = transpose(smooth(PK_profile(1:5:end)));

    figure(1)
    subplot('position',[0.06 0.55 0.75 0.4])
    %subplot('position',[left bottom width height])
    %plot(x/1000,zb,'k','linewidth',2)
    %hold on
    plot(x/1000,zb_ccp+H,'c','linewidth',3) % updates for topography
    hold on
    plot(x/1000,zb_ccp,'k','linewidth',2)
    legend('Glacier','Pine Creek Topography after smoothing operation')
    %plot(x/1000,z,'linewidth',1) % normal glacier on black topo
    plot(x/1000,ELA0*ones(size(x)),'g--','linewidth',2.5)
    plot(x/1000,(ELA0+290)*ones(size(x)),'g--','linewidth',1)
    plot(x/1000,(ELA0-290)*ones(size(x)),'g--','linewidth',1)
    axis([0.2 30 2500 4200])
    title('Pine Creek Valley paleoglacier, LGM numerical reconstruction') 
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('elevation [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    hold off
    
    %plot(x/1000,zb_ccp+H,'-','linewidth',0.25) % updates for topography

    
    subplot('position',[0.88 0.55 0.1 0.4])
            %subplot('position',[left bottom width height])
    plot(b,z,'k','linewidth',2)
    hold on
    plot(b0,zb,'k--','linewidth',2.5)
    plot(zeros(size(z)),z,'g','linewidth',2.5)
    axis([minb 1.5*bcap 2500 4200])
    title('mass balance')
    xlabel('b(z) [m/yr]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
  
    time=num2str(t(i));
    timetext=strcat(time,' yrs')
    text(-42,4000,timetext,'fontsize',18)
    hold off
    
    Qanal = cumsum(b)*dx; % analytic solution for ss ice discharge
    Qanal = max(Qanal,0);
    
    %figure(3)
    subplot('position',[0.1 0.1 0.38 0.3])
 %  subplot('position',[left bottom width height])
    plot(xedge/1000,Q/1000,'c','linewidth',3)
    hold on
    plot(x/1000,Qanal/1000,'k--','linewidth',1.5)
    axis([0 30 0 30])
    legend('Discharge','Steady state analytic solution')
    title('steady state ice discharge')
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('ice discharge [10^3 m^2/yr]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    hold off
    
    %figure(4)
    subplot('position',[0.58 0.1 0.38 0.3])
 %  subplot('position',[left bottom width height])
    plot(x/1000,H,'linewidth',0.5)
    hold on
    axis([0 30 0 450])
    title('cumulative glacier ice thickness')
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('ice thickness [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial') 
    pause(0.01)
    
end

end
%% finalize

figure(2)
subplot(2,2,1)
plot(t/1000,ELA,'g','linewidth',0.25)
hold on
plot(t/1000,(ELA0)*ones(size(t)),'g--','linewidth',2.5)
plot(t/1000,(ELA0+290)*ones(size(t)),'g--','linewidth',1.5)
plot(t/1000,(ELA0-290)*ones(size(t)),'g--','linewidth',1.5)
    xlabel('time [ka]','fontname','arial','fontsize',18)
    ylabel('ELA [m]','fontname','arial','fontsize',18)
    set(gca,'XDIR','reverse','fontsize',18,'fontname','arial')
    axis([0 4 ELA0-300 ELA0+600])
subplot(2,2,2)
    plot(t/1000,dam_ht,'c','linewidth',3)
    xlabel('time [ka]','fontname','arial','fontsize',18)
    ylabel('dam height [m]','fontname','arial','fontsize',18)
    set(gca,'XDIR','reverse','fontsize',18,'fontname','arial')

load time_10_26.txt
load ELA_normalized.txt
    
x1 = transpose(time_10_26(:,1));
y1 = transpose(ELA_normalized(:,1));

dxx = 0.01;
dyy = 0.01;

xx = 0:dxx:16;
yy = 0:dyy:16;

t_1 = interp1(x1,xx,'spline');
z_1 = interp1(y1,yy,'spline');
    
    subplot(2,2,4)
    plot(t_1/1000,z_1+100,'g','linewidth',3)
    axis([12 26 2800 5600])
    hold on
plot(t_1/1000,(ELA0)*ones(size(t_1)),'g--','linewidth',2.5)
plot(t_1/1000,(ELA0+290)*ones(size(t_1)),'g--','linewidth',1.5)
plot(t_1/1000,(ELA0-290)*ones(size(t_1)),'g--','linewidth',1.5)
    xlabel('time [ka]','fontname','arial','fontsize',18)
    ylabel('ELA [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    
    subplot(2,2,3)
    plot(t_1/1000,z_1+100,'g','linewidth',3)
    hold on
plot(t_1/1000,(ELA0)*ones(size(t_1)),'g--','linewidth',2.5)
plot(t_1/1000,(ELA0+290)*ones(size(t_1)),'g--','linewidth',1.5)
plot(t_1/1000,(ELA0-290)*ones(size(t_1)),'g--','linewidth',1.5)
    axis([17 21 ELA0-300 ELA0+600])
    xlabel('time [ka]','fontname','arial','fontsize',18)
    ylabel('ELA [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')

%% end


