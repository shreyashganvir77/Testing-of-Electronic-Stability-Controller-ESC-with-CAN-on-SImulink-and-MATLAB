clear;
clc;

% reading steering angle and velocity parameter
load('Drv_DeltaSteer.mat');
load('Veh_Vx.mat');
drivingProfile = csvread('parameter.csv',1,0);

% extracting the sampling instances
n = floor(drivingProfile(9)/0.005);
time = floor(size(Drv_DeltaSteer,1)/n);
delete 'Time.csv'
fileTime = fopen('Time.csv','w+');
fprintf(fileTime,'%10.20f\n',time);
fclose(fileTime);

steeringAngle = zeros(1, time);
velocity = zeros(1, time);
length = size(Drv_DeltaSteer,1);
j = 1;
for i=1:length
    if mod(i, n)==0
        steeringAngle(j) = Drv_DeltaSteer(i);
        velocity(j) = Veh_Vx(i);
        j = j + 1;
    end
end

% parameters
m = drivingProfile(1); %kg
v = max(Veh_Vx); % m/s
lf = drivingProfile(2); %m
lr = drivingProfile(3); %m
l = lr + lf;
Cf = drivingProfile(4); %N/rad
Cr = drivingProfile(5); %N/rad
Iz = drivingProfile(6); %kg m^2

mu = drivingProfile(7);
g = 9.8; %m/s^2
steering_ratio = drivingProfile(8); %15:1
Ku = (m*((lr*Cr)-(lf*Cf)))/(l*Cf*Cr);
max_sideslip = atan(0.02*mu*g);

a11 = (-2)*((Cf+Cr)/(m*v));
a12 = -1-(((2*Cf*lf)-(2*Cr*lr))/(m*v*v));
a21 = (-2)*(((lf*Cf)-(lr*Cr))/Iz);
a22 = (-2)*(((lf*lf*Cf)+(lr*lr*Cr))/(Iz*v));
b1 = (2*Cf)/(m*v);
b2 = (2*Cf*lf)/Iz;
b3 = 0;
b4 = 1/Iz;

% Discretizing the state space
A = [a11 a12;
     a21 a22];
B = [b1 b3;
     b2 b4];
C = [0 1;
    (v*a11) (v*(a12 + 1))];
D = [0 0;
    (v*b1) 0];

sys = ss(A,B,C,D);
Ts = drivingProfile(9); % Sampling period
sys_d = c2d(sys,Ts, 'zoh');

A = sys_d.A;
B = sys_d.B;
C = sys_d.C;
D = sys_d.D;

% Feedback gain computatioon
co = ctrb(sys_d);
isControllable = [rank(A) == rank(co)]

K = [   5.261, -0.023;
     -414911.257,   57009.478];
 
% Observer gain computation
ob = obsv(sys_d);
isObservable = [rank(A) == rank(ob)]
QN = 1;
RN = eye(size(C));
[KEST,L,P,M,Z] = kalman(sys_d,QN,RN);
abs(eig(A - L*C))

% Writing K, L values to csv files
delete 'K.csv'
delete 'L.csv'
fileK = fopen('K.csv','w+');
fileL = fopen('L.csv','w+');
fprintf(fileK,'%10.20f %10.20f\n',K.');
fprintf(fileL,'%10.20f %10.20f\n',L.');
fclose(fileK);
fclose(fileL);

plot_sideslip = zeros(1, time);

plot_yaw = zeros(1, time);
plot_sidesliphat = zeros(1, time);
plot_yawhat = zeros(1, time);

x0 = [0;0]; 
x = x0;
xhat0 = [0;0];
xhat = xhat0;
y = C*x;
yhat = C*xhat;
u0 = [0;0]; 
u = u0;

delete 'A.csv'
delete 'B.csv'
delete 'C.csv'
delete 'D.csv'
delete 'X.csv'
delete 'Z.csv'
delete 'Y.csv'
delete 'U.csv'
delete 'LatDist.csv'
delete 'Velocity.csv'
delete 'Steering.csv'
fileA = fopen('A.csv','a+');
fileB = fopen('B.csv','a+');    
fileC = fopen('C.csv','a+');
fileD = fopen('D.csv','a+');
fileVelocity = fopen('Velocity.csv','a+');
fileSteering = fopen('Steering.csv','a+');
fileX = fopen('X.csv','a+');
fileXHat = fopen('Z.csv','a+');
fileY = fopen('Y.csv','a+');
fileU = fopen('U.csv','a+');

fprintf(fileX,'%10.10f %10.10f\n',x.');
fprintf(fileXHat,'%10.10f %10.10f\n',xhat.');
fprintf(fileY,'%10.10f %10.10f\n',y.');
fprintf(fileU,'%10.10f %10.10f\n',u.');

for i=1:time
    delta = steeringAngle(i)*(pi/180)*steering_ratio;
    v = velocity(i);

    fprintf(fileSteering,'%10.10f\n',delta);
    fprintf(fileVelocity,'%10.10f\n',v);    
    
    a11 = (-2)*((Cf+Cr)/(m*v));
    a12 = -1-(((2*Cf*lf)-(2*Cr*lr))/(m*v*v));
    a21 = (-2)*(((lf*Cf)-(lr*Cr))/Iz);
    a22 = (-2)*(((lf*lf*Cf)+(lr*lr*Cr))/(Iz*v));
    b1 = (2*Cf)/(m*v);
    b2 = (2*Cf*lf)/Iz;
    b3 = 0;
    b4 = 1/Iz;

    A = [a11 a12;
         a21 a22];
    B = [b1 b3;
         b2 b4];
    C = [0 1;
        (v*a11) (v*(a12 + 1))];
    D = [0 0;
        (v*b1) 0];
    
    sys = ss(A,B,C,D);
    Ts = 0.04;
    sys_d = c2d(sys,Ts, 'zoh');
    
    A = sys_d.A;    
    fprintf(fileA,'%10.10f %10.10f\n',A.');
    B = sys_d.B;
    fprintf(fileB,'%10.10f %10.10f\n',B.');
    C = sys_d.C;
    fprintf(fileC,'%10.10f %10.10f\n',C.');
    D = sys_d.D;
    fprintf(fileD,'%10.10f %10.10f\n',D.');
    
    a = A(1,2)/A(2,2);
    b = B(1,1) - ((B(2,1)*A(1,2))/A(2,2));
    c = v/(A(1,1) - ((A(2,1)*A(1,2))/A(2,2)));    
    
    r = y - yhat;
    x = A*x + B*u;
    xhat = A*xhat + B*u + L*r;
    u = -(K*xhat) + [delta;0];
    y = C*x + D*u;
    yhat = C*xhat + D*u;
    
    fprintf(fileX,'%10.10f %10.10f\n',x.');
    fprintf(fileXHat,'%10.10f %10.10f\n',xhat.');
    fprintf(fileY,'%10.10f %10.10f\n',y.');
    fprintf(fileU,'%10.10f %10.10f\n',u.');
    
    plot_sideslip(i) = x(1);
    plot_sidesliphat(i) = xhat(1);
    plot_yaw(i) = x(2);
    plot_yawhat(i) = xhat(2);
end

fontsize = 10;
linewidth = 1;

clf;
subplot(1,2,1);
hold on;
plot(plot_yaw,'Linewidth',linewidth);
plot(plot_yawhat,'Linewidth',linewidth);
set(gca,'FontSize',fontsize)
xlabel('Time(x40x10^{-3})(s)','FontSize',fontsize);
ylabel('rad/s','Fontsize',fontsize);
title('yaw rate');
grid on;
hold off;

subplot(1,2,2);
hold on;
plot(plot_sideslip,'Linewidth',linewidth);
plot(plot_sidesliphat,'Linewidth',linewidth);
set(gca,'FontSize',fontsize)
xlabel('Time(x40x10^{-3})(s)','FontSize',fontsize);
ylabel('rad','Fontsize',fontsize);
title('sideslip');
grid on;
hold off;

fclose(fileA);
fclose(fileB);
fclose(fileC);
fclose(fileD);
fclose(fileSteering);
fclose(fileVelocity);
fclose(fileX);
fclose(fileXHat);
fclose(fileY);
fclose(fileU);

% moving preprocessed files to folder 'preprocess'
mkdir preprocess
movefile Time.csv preprocess
movefile K.csv preprocess
movefile L.csv preprocess
movefile A.csv preprocess
movefile B.csv preprocess
movefile C.csv preprocess
movefile D.csv preprocess
movefile X.csv preprocess
movefile Z.csv preprocess
movefile Y.csv preprocess
movefile U.csv preprocess
movefile Velocity.csv preprocess
movefile Steering.csv preprocess