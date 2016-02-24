%%%  DATE: Feb. 23rd 2016
%%   CLEANED UP GLACIER MODEL FOR CLEAR CREEK - STAMPING GLACIAL LINES
%    1D FTCS STAGGERED GRID NUMERICAL MODEL: CLEAR CREEK GLACIER – SI units

%    COLE C. PAZAR    and    ROBERT S. ANDERSON

    clearvars

    figure(1)  % main figure for animation of the glacier, full screen !
    clf

    figure(2)  % secondary climate figure after animation has ended
    clf


%% initialize

%  constants

step  = 200;        % matrix sizes for the model
rho_i = 917;        % density of glacial ice
g     = 9.81;       % gravitational acceleration near the surface
A     = 2.1e-16;    % glenn-nye flow law parameter [=] Pa-3 yr-1
slide = 0.5;        % ratio of sliding speed to internal deformation speed

%  set up distance array

    xmax  = 30000;
    
       dx = xmax/step;

        x = dx/2:dx:xmax-(dx/2);

    xedge = 0:dx:xmax;

%  now for the Clear Creek, Colorado case

    zb0    = 3425;
    slope0 = 0.0285;       
    zb1    = 516;
    xstar  = 900;           
    zbline = zb0 - (slope0*x);
    zbexp  = zb1 * exp(-x/xstar);

    zb = zbline + zbexp;

%       valley width as a function of distance downvalley
%       turn this off by setting Wmin=W0;

    W0    = 3000;
    Wstar = 4000;
    Wmin  = 1700;
    W     = Wmin + (W0-Wmin)*exp(-x./Wstar);

 Wedge = W(1:end-1)+0.5*diff(W); % interpolates valley width to cell edges

    Wedge = [Wedge(1) Wedge Wedge(end)];

%   now put in a damming wall at the bottom of the valley

    slopewall = 0.3;
    xwall     = 29450;
    zbwall    = 2800 + (slopewall * (x-xwall));
    
        zb = max(zb,zbwall);
        
    beyond     = find(x>29450);
    slope2     = 0.01;
    zb(beyond) = 2900+(slope2*(x(beyond)-xwall));

    dambase  = find(zb==min(zb));
    xdambase = x(dambase);
    zbmin    = min(zb);
    zbmax    = max(zb);

    H = zeros(size(x));
    z = zb+H;

%  meteorology and mass balance

    ELA0 = 3300;      % SET THE AVERAGE ELA ––> SEE BRUGGER (2010)
    
    sigma_ELA = 300;  % uncertainty in the ELA

    dbdz = 0.006;          % m/y/m, typically 0.01
    bcap = 1.25;           % m/yr, usually 1.25-2.00
    b0 = dbdz*(z-ELA0);
    b0 = min(b0,bcap);
    minzb = find(zb==min(zb));
    minb = b0(minzb);

%  set up time array

    dt   = 0.0025; % time step has to be small for glaciers
    tmax = 4000;   % max time interval of growth, 
    
        t = 0:dt:tmax;

%  plotting controls 

    imax = length(t);
    nplots = 80;
    tplot = tmax/nplots;
    nframe = 0;

%  new way to control the Climate:

    big_period = 100000; % years for the 100 ka milankovich cycle
    
    med_period = 2500;   % years so that the peak is centered

    sol_period = 11;     % years for solar cycle

    big_shift = -50000;  % shift the periods

    med_shift = 1000;    % shift the periods

sol_shift = 11*randn(size(t));  % include a solar cycle  w/  1 m / yr

randomsize_t = randn(size(t));  % randomize variables

    big = (sigma_ELA/4)*sin(2*pi*(t+big_shift)/big_period);
    medium = (sigma_ELA/4).*sin(2*pi*(t+med_shift)/med_period);
    solar = (sol_shift).*sin(2*pi*t/sol_period);

ELA = ELA0 + sigma_ELA + big + medium + solar - 180;  % final climate


%% run the model

for i = 1:imax
    
b = dbdz*(z-ELA(i)); % local net balance calculated at cell centers at ice surface
b = min(b,bcap);

Hedge = H(1:end-1)+0.5*diff(H); % interpolates ice thickness to cell edges
S = abs(diff(z)/dx); % slope of ice surface calculated at cell edges

Udef = (A/5).*((rho_i*g*S).^3).*(Hedge.^4); %mean defm speed
Q = (A/5).*((rho_i*g*S).^3).*(Hedge.^5); % ice discharge due to internal deformation
Qsl = slide * Udef.*Hedge;
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
        
    load CC_new_profile.txt
    zb_ccp = transpose(smooth(CC_new_profile(1:5:end)));

    figure(1)
    subplot('position',[0.06 0.55 0.75 0.4])
    plot(x/1000,zb_ccp+H,'linewidth',1) % updates the glacier for topography
    hold on
    plot(x/1000,zb_ccp,'k','linewidth',2)
    legend('Glacier','Clear Creek Topography after smoothing operation')
    plot(x/1000,ELA0*ones(size(x)),'g--','linewidth',2.5)
    plot(x/1000,(ELA0+290)*ones(size(x)),'g--','linewidth',1)
    plot(x/1000,(ELA0-290)*ones(size(x)),'g--','linewidth',1)
    axis([0.2 30 2500 4200])
    title('Clear Creek Valley paleoglacier, LGM numerical reconstruction') 
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('elevation [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    
    subplot('position',[0.88 0.55 0.1 0.4])
    plot(b,z,'k','linewidth',2)
    hold on
    plot(b0,zb,'k--','linewidth',2.5)
    plot(zeros(size(z)),z,'g','linewidth',2.5)
    axis([minb 1.5*bcap 2500 4200])
    xlabel('b(z) [m/yr]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
  
    time=num2str(t(i));
    timetext=strcat(time,' yrs')
    text(-42,4000,timetext,'fontsize',18)
    hold off
    
    Qanal = cumsum(b)*dx; % analytic solution for ss ice discharge
    Qanal = max(Qanal,0);
    
    subplot('position',[0.1 0.1 0.38 0.32])

    plot(xedge/1000,Q/1000,'linewidth',1)
    hold on
    plot(x/1000,Qanal/1000,'k--','linewidth',1.5)
    axis([0 30 0 30])
    legend('Discharge','Steady state analytic solution')
    title('Steady state ice discharge')
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('ice discharge [10^3 m^2/yr]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    
    subplot('position',[0.58 0.1 0.38 0.32])
    plot(x/1000,H,'linewidth',0.5)
    hold on
    axis([0 30 0 450])
    title('Cumulative glacier ice thickness')
    xlabel('horizontal distance [km]','fontname','arial','fontsize',18)
    ylabel('ice thickness [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial') 
    pause(0.01)
    
end

end
%% finalize

    figure(2)
    
    subplot(2,2,1)
    title('Climate reconstruction')
    plot(t/1000,ELA,'g','linewidth',0.25)
    hold on
    plot(t/1000,(ELA0)*ones(size(t)),'g--','linewidth',2.5)
    plot(t/1000,(ELA0+290)*ones(size(t)),'g--','linewidth',1.5)
    plot(t/1000,(ELA0-290)*ones(size(t)),'g--','linewidth',1.5)
    xlabel('time [ka]','fontname','arial','fontsize',18)
    ylabel('ELA [m]','fontname','arial','fontsize',18)
    set(gca,'XDIR','reverse','fontsize',18,'fontname','arial')
    axis([0 4 ELA0-300 ELA0+600])
    
    subplot(2,2,3)
    title('Tracking the ice dam')
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

t_1 = interp1(x1,xx,'spline'); % splines the data as an interpolation for
z_1 = interp1(y1,yy,'spline'); % fitting the ELA curve onto the plots
    
    subplot(2,2,2)
    title('Scaled ∆O^18 curve')
    plot(t_1/1000,z_1+150,'g','linewidth',3)
    hold on
    plot(t_1/1000,(ELA0)*ones(size(t_1)),'g--','linewidth',2.5)
    plot(t_1/1000,(ELA0+290)*ones(size(t_1)),'g--','linewidth',1.5)
    plot(t_1/1000,(ELA0-290)*ones(size(t_1)),'g--','linewidth',1.5)
    axis([17 21 ELA0-300 ELA0+600])
    xlabel('time [ka]','fontname','arial','fontsize',18)
    ylabel('ELA [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')
    
    
    subplot(2,2,4)
    title('Full scaled ∆O^18 curve')
    plot(t_1/1000,z_1+150,'g','linewidth',3)
    axis([12 26 2800 5600])
    hold on
    plot(t_1/1000,(ELA0)*ones(size(t_1)),'g--','linewidth',2.5)
    plot(t_1/1000,(ELA0+290)*ones(size(t_1)),'g--','linewidth',1.5)
    plot(t_1/1000,(ELA0-290)*ones(size(t_1)),'g--','linewidth',1.5)
    xlabel('time [ka]','fontname','arial','fontsize',18)
    ylabel('ELA [m]','fontname','arial','fontsize',18)
    set(gca,'fontsize',18,'fontname','arial')


%% end

        %  subplot('position',[left bottom width height])
