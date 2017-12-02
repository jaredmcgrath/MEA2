g = 9.81; %Acceleration by gravity
dt = 0.025; %time step increment
vWind = [4 0]; %wind speed vector CHANGE THIS VALUE TO [10 0] or [-10 0]
mass = 4; %mass of projectile
rho = 1.225; %density of air
coeffDrag = 0.5; %Coefficient of drag, generally 0.5 for a sphere
area = 0.097986; %cross-sectional area of projectile
maxDist = 0; %var for max distance
maxTheta = 0; %var for max theta
maxV = 0;

%define trebuchet constants
mBlock = 133*mass; %mass of block, kg
I = 4453.93; %moment of intertia for trebuchet
lString = 7; %length of string, meters
lRod = 7; %length of rod, meters
hBlock = 11.39; %height of block, meters


%find optimal theta
for theta = 0:0.05:90
    a = sqrt((2 * (lString + lRod)^2) - (2 * (lRod + lString)^2 * cosd(theta)));
    hp = a * sind(45 + asind( ((lRod + lString)*sind(theta))/a ));
    x = -(lString + lRod) * cosd(theta - 45);
    vi = sqrt( (2*((mBlock * g * hBlock) + (mass * g * hp)))/(mass + I/((lRod + lString)^2)) );
    v=[vi*cosd(theta) vi*sind(theta)]; %velocity into components
    p = [x hp]; %update this to be the initial position vector
    while p(2)>0
        vLast = v;

        %      1/2 *density  *  Cd  *  x-sect-area * v^2   *   opposite sign of current velocity 
        Fdrag = 0.5 * rho * coeffDrag * area * (vLast.^2) .* (-1 * sign(vLast));

        %Fnet is gravity plus drag
        Fnet = [0 -g*mass] + Fdrag;

        a = Fnet / mass; %acceleration
        v = vLast + a*dt; %new current velocity

        vAvg = 0.5 * (vLast + v); %avg velocity over last time interval
        p = p + (vAvg+vWind)*dt; %update position
    end
    if p(1)>maxDist %check to find max distance/theta
        maxDist = p(1);
        maxTheta = theta;
        maxV = vi;
    end
    disp(['Distance travelled: ' num2str(p(1))]);
end
disp(['Max v: ' num2str(maxV)]);
disp(['Max distance: ' num2str(maxDist)]);
disp(['Angle: ' num2str(maxTheta)]);

%draw plot for best theta
figure
theta=maxTheta;
v=[vi*cosd(theta) vi*sind(theta)];
p = [0 2];
while p(2)>0
    vLast = v;

    Fdrag = 0.5 * rho * coeffDrag * area * (vLast.^2) .* (-1 * sign(vLast));

    Fnet = [0 -g*mass] + Fdrag;

    a = Fnet / mass;
    v = vLast + a*dt;

    vAvg = 0.5 * (vLast + v);
    p = p + (vAvg+vWind)*dt;

    hold on
    plot(p(1),p(2),'r.');
    axis equal
    getframe;
end