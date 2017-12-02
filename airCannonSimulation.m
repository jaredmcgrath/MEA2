%define constants. All quantities are SI units
r = 0.12488; %radius, m
Vp = 1/3 * pi * r^3; %constant volume under pumpkin, m^3
T = 293.15; %temperature, K
Vres = 10; %volume of reservoir, m^3
Pext = 101300; %external pressure, Pa
g = 9.81; %gravity, N/kg
m = 4; %mass of pumpkin, kg
Asurf = 2 * pi * r^2; %half surface area of pumpkin, m^2
L = 12; %length of barrel
coeffFlow = 87/100; %flow rate coefficient, m^3/(s*Pa)
R = 8.314; %gas constant, Pa*m^3/(mol*K)
maxOutlet = 1e+7; %max outlet flow, Pa
Pmax = 2500000; %max pressure on pumpkin, Pa. Based on estimates that are well-established in the field of punkin chunkin
time = 0;

%Copy-pasting program indSimulation.m
aGrav = 9.81; %Acceleration by gravity
vWind = [4 0]; %wind speed vector CHANGE THIS VALUE EMMA. should be [10 0] or [-10 0]
mass = 4; %mass of projectile
rho = 1.225; %density of air
coeffDrag = 0.5; %Coefficient of drag, generally 0.5 for a sphere
area = 0.097986; %cross-sectional area of projectile
maxDist = 0; %var for max distance
maxTheta = 0; %var for max theta

%data point tracking variables
times = [];
pData = [];
dData = [];
vData = [];
aData = [];

%find optimal theta
for theta = 0:0.5:90
    %define dynamic variables for pressure sim
    dt = 0.0001; %change in time for sim, s
    Pbar = 101300; %initial pressure in barrel, Pa
    Vbar = Vp; %volume under pumpkin in barrel, m^3
    nBar = (Pbar*Vbar)/(R*T); %number of moles in barrel initially
    d = 0; %displacement of pumpkin in barrel, m
    v = 0; %velocity of pumpkin, m/s
    a = 0; %acceleration of pumpkin, m/s^2
    Pout = 1; %pressure of outflow stream, Pa
    time = 0;
    
    while d<L
        vPrev = v;

        PoutPrev = Pout;
        nBarPrev = nBar;
        PbarPrev = Pbar;
        VbarPrev = Vbar;

        %calculate small change in gas values
        dn = Pout*coeffFlow*dt/(R*T);

        %update gas values in barrel
        Vbar = Vp + (pi * r^2 * d); %new volume in barrel
        nBar = nBarPrev + dn;
        Pbar = nBar*R*T/Vbar;

        %take average values between simulation points
        VbarAvg = (VbarPrev + Vbar)/2; %avg volume between simulation points
        PbarAvg = (PbarPrev + Pbar)/2;
        nBarAvg = (nBarPrev + nBar)/2;

        Pout = Pmax-(Pbar + Pout*coeffFlow*dt/VbarAvg);

        if Pout > maxOutlet
            Pout = maxOutlet;
        end


        Fnet = Asurf * (PbarAvg - Pext) - m*g*sind(theta);
        a = Fnet/m;
        v = vPrev + a*dt;
        vAvg = (v + vPrev)/2;
        d = d + vAvg * dt;
        if d<0
            d=0;
        end

        time = time + dt;
    end

    disp(['Exit velocity: ' num2str(v)]);
    
    %at this point, exit velocity has been found, passing to wind sim part
    %wind sim part
    dt = 0.025; %time step increment
    vv=[v*cosd(theta) v*sind(theta)]; %velocity into components
    p = [0 L*sind(theta)]; %update this to be the initial position vector
    while p(2)>0
        vLast = vv;

        %      1/2 *density  *  Cd  *  x-sect-area * v^2   *   opposite sign of current velocity 
        Fdrag = 0.5 * rho * coeffDrag * area * (vLast.^2) .* (-1 * sign(vLast));

        %Fnet is gravity plus drag
        Fnet = [0 -aGrav*mass] + Fdrag;

        a = Fnet / mass; %acceleration
        vv = vLast + a*dt; %new current velocity

        vAvg = 0.5 * (vLast + vv); %avg velocity over last time interval
        p = p + (vAvg+vWind)*dt; %update position
    end
    if p(1)>maxDist %check to find max distance/theta
        maxDist = p(1);
        maxTheta = theta;
        vExit = v;
    end
    disp(['Distance travelled: ' num2str(p(1))]);
end
disp(['Exit Velocity final: ' num2str(vExit)]);
disp(['Max distance: ' num2str(maxDist)]);
disp(['Angle: ' num2str(maxTheta)]);

%draw plot for best theta, both wind sim and pressure sim
%pressure sim

dt = 0.0001;

%define dynamic variables for pressure sim
Pres = 15000000; %initial pressure in reservoir, Pa
Pbar = 101300; %initial pressure in barrel, Pa
Vbar = Vp; %volume under pumpkin in barrel, m^3
nBar = (Pbar*Vbar)/(R*T); %number of moles in barrel initially
d = 0; %displacement of pumpkin in barrel, m
v = 0; %velocity of pumpkin, m/s
a = 0; %acceleration of pumpkin, m/s^2
Pout = maxOutlet; %pressure of outflow stream, Pa
theta=maxTheta;
vAvg = 0; %avergae velocity
moved = false; %has the pumpkin moved yet?
tMove = 0;
pMove = 0;
time=0;

%configure figures
figure(1)
xlabel('Time, s'), ylabel('Pressure, Pa')
title('Pressure in Cannon Barrel over Time')
figure(2)
xlabel('Time, s'), ylabel('Acceleration, m/s^2')
title('Acceleration of Pumpkin over Time')
figure(3)
xlabel('Time, s'), ylabel('Velocity, m/s')
title('Velocity of Pumpkin over Time')
figure(4)
xlabel('Time, s'), ylabel('Displacement, m')
title('Displacement of Pumpkin over Time')

while d<L
    vPrev = v;
    
    %record data points
    pData(end+1) = Pbar;
    aData(end+1) = a;
    vData(end+1) = vAvg;
    dData(end+1) = d;
    times(end+1) = time;
    
    PoutPrev = Pout;
    nBarPrev = nBar;
    PbarPrev = Pbar;
    VbarPrev = Vbar;
    
    %calculate small change in gas values
    dn = Pout*coeffFlow*dt/(R*T);
    
    %update gas values in barrel
    Vbar = Vp + (pi * r^2 * d); %new volume in barrel
    nBar = nBarPrev + dn;
    Pbar = nBar*R*T/Vbar;
    
    %take average values between simulation points
    VbarAvg = (VbarPrev + Vbar)/2; %avg volume between simulation points
    PbarAvg = (PbarPrev + Pbar)/2;
    nBarAvg = (nBarPrev + nBar)/2;
    
    Pout = Pmax-(Pbar + Pout*coeffFlow*dt/VbarAvg);
    
    if Pout > maxOutlet
        Pout = maxOutlet;
    end
    if d>0 &&  ~moved
        tMove = time;
        pMove = Pbar;
        moved = true;
    end
    
    Fnet = Asurf * (PbarAvg - Pext) - m*g*sind(theta);
    a = Fnet/m;
    v = vPrev + a*dt;
    vAvg = (v + vPrev)/2;
    d = d + vAvg * dt;
    if d<0
        d=0;
    end
    

    time = time + dt;
end
figure(1)
hold on
plot(times',pData');
txt1 = '\leftarrow Pumpkin Begins to Move';
text(tMove,pMove,txt1);
figure(2)
hold on
plot(times',aData');
figure(3)
hold on
plot(times',vData');
figure(4)
hold on
plot(times',dData');
getframe;
disp(['Exit velocity: ' num2str(v)]);
disp(['time: ' num2str(time)]);


%wind sim part

dt = 0.025; %time step increment
figure(5)
xlabel('Horizontal Distance, m'), ylabel('Vertical Height, m')
title('Position of Pumpkin After Launch over Time')

v=[v*cosd(theta) v*sind(theta)];
p = [0 L*sind(theta)];
while p(2)>0
    vLast = v;

    Fdrag = 0.5 * rho * coeffDrag * area * (vLast.^2) .* (-1 * sign(vLast));

    Fnet = [0 -aGrav*mass] + Fdrag;

    a = Fnet / mass;
    v = vLast + a*dt;

    vAvg = 0.5 * (vLast + v);
    p = p + (vAvg+vWind)*dt;

    figure(5)
    hold on
    plot(p(1),p(2),'r.');
    axis equal
end
getframe;