% 1D FTCS STAGGERED GRID NUMERICAL MODEL: CLEAR CREEK GLACIER
% COLE C. PAZAR AND ROBERT S. ANDERSON
clearvars
figure(1)
clf
% figure(2)
% clf
figure(3)
clf
% figure(4)
% clf

%% initialize
step = 500;
rho_i = 917;
g = 9.81;
A = 2.1e-16; % Pa-3 yr-1
slide_ratio = 0.5;% ratio of sliding speed to internal defm speed

% initial topography
% zbmax = 4000;
% zbmin = 2000;
% xstar = 7000;

xmax = 28000;
dx = xmax/step;
x = dx/2:dx:xmax-(dx/2);
xedge = 0:dx:xmax;

%zb = zbmin+(zbmax-zbmin)*exp(-x/xstar);

% now a clear creek colorado case
zb0 = 3500;
slope0 = 0.0338;
%slope0 = 0.03;
zb1 = 516;
xstar = 1687;
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
slopewall = 0.4;
xwall = 26000;
zbwall = 2680 + (slopewall * (x-xwall));
zb = max(zb,zbwall);
beyond = find(x>26500);
slope2 = 0.1;
zb(beyond) = 2800+(slope2*(x(beyond)-25500));

%     load clear_creek_profile.txt
% 
%     zb = transpose(flipud(clear_creek_profile(1:2:end,1:2:end,1:2:end)));

dambase = find(zb==min(zb));
xdambase = x(dambase);
zbmin = min(zb);
zbmax = max(zb);

H = zeros(size(x));
%Hedge = zeros(size(xedge));
z = zb+H;

% now set up meteorology (some based upon brugger's suggestion in sawatch
% paper 2006)
%ELA0 = 3300;
ELA0 = 3270;  % AVERAGE: SEE BRUGGER (2010)

sigma_ELA = 290;  % used to be 250

dbdz = 0.006; % m/y/m   0.01
bcap = 1.5; % m/yr (usually 2.0 m/yr or 1.25)
b0 = dbdz*(z-ELA0);
b0 = min(b0,bcap);
minzb = find(zb==min(zb));
minb = b0(minzb);

dt = 0.0025;
tmax = 1000;

t = 0:dt:tmax;
imax = length(t);
nplots = 20;
tplot = tmax/nplots;

period = 100; % years
ELA = ELA0 + (sigma_ELA)*randn(size(t)); % + (sigma_ELA/2)*sin(2*pi*t/period);
%ELA = ELA0 + sigma_ELA*sin(2*pi*t/period);
%ELA = 3300*ones(size(t));

nframe = 0;

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
term(i) = x(glacier(end));

% track dam height
if(x(glacier(end))>xdambase)
    dam_ht(i) = zb(glacier(end))-min(zb);
else
    dam_ht(i) = 0;
end
    dam_ht(i) = max(0,dam_ht(i));
    
% now for some plotting
if rem(t(i),tplot)==0
    nframe = nframe + 1
    %a = 0.971;% scaling factor 1
    %load clear_creek_profile.txt
    %zb_ccp = transpose(flipud(clear_creek_profile(1:2:end,1:10:end,1:10:end)));
    figure(1)
    %subplot(1,3,1)
    subplot('position',[0.06 0.55 0.75 0.4])
        %subplot('position',[left bottom width height])
    plot(x/1000,zb,'k','linewidth',2)
    % legend('approximated bed topography')
    hold on
    % plot(x/1000,zb_ccp./a,'r','linewidth',1.5)
    % fill(x/1000,H,'c')
    plot(x/1000,z,'linewidth',1)
    plot(x/1000,ELA0*ones(size(x)),'g--','linewidth',2.5)
    plot(x/1000,(ELA0+290)*ones(size(x)),'g--','linewidth',1)
    plot(x/1000,(ELA0-290)*ones(size(x)),'g--','linewidth',1) % analytic envelopes of precipitation window
    axis([0 28 2600 4200])
    title('Clear Creek Valley paleoglacier, LGM numerical reconstruction: ELA probability = 3270 ± 290 m') 
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('elevation [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    
    %subplot(1,3,2)
    subplot('position',[0.88 0.55 0.1 0.4])
            %subplot('position',[left bottom width height])
    plot(b,z,'k','linewidth',2)
    hold on
    plot(b0,zb,'k--','linewidth',2.5)
    plot(zeros(size(z)),z,'g','linewidth',2.5)
    axis([minb 1.5*bcap 2600 4200])
    title('mass balance')
    xlabel('b(z) [m/yr]','fontname','arial','fontsize',18)
    %ylabel('elevation (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    % now stamp times
    time=num2str(t(i));
    timetext=strcat(time,' yrs')
    text(-3,4000,timetext,'fontsize',18)
    hold off
    
    Qanal = cumsum(b)*dx; %analytic solution for ss ice discharge
    Qanal = max(Qanal,0);
    
    %figure(3)
    subplot('position',[0.1 0.1 0.38 0.3])
 %subplot('position',[left bottom width height])
    plot(xedge/1000,Q/1000,'linewidth',1.5)
    hold on
    plot(x/1000,Qanal/1000,'k--','linewidth',1)
    axis([0 30 0 30])
    title('steady state ice discharge')
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('ice discharge [10^3 m^2/yr]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')    
    %figure(4)
    subplot('position',[0.58 0.1 0.38 0.3])
 %subplot('position',[left bottom width height])
    plot(x/1000,H)
    hold on
    axis([0 30 0 450])
    title('cumulative ice thickness')
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('ice thickness [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    pause(0.1)
end

end
%% finalize
figure(3)
% subplot(3,1,1)
% plot(t,term/1000)
%     xlabel('time [years]','fontname','arial','fontsize',18)
%     ylabel('terminus position [km]','fontname','arial','fontsize',18)
%     set(gca,'fontsize',18,'fontname','arial')
% subplot(3,1,2)
% plot(t,dam_ht)
%     xlabel('time [years]','fontname','arial','fontsize',18)
%     ylabel('dam height [m]','fontname','arial','fontsize',18)
%     set(gca,'fontsize',18,'fontname','arial')
% subplot(3,1,3)
plot(t,ELA)
hold on
plot(t,(ELA0+290)*ones(size(x)),'g--','linewidth',1)
plot(t,(ELA0-290)*ones(size(x)),'g--','linewidth',1) % analytic envelopes of precipitation window
    xlabel('time [years]','fontname','arial','fontsize',18)
    ylabel('ELA [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    
%% end of code
