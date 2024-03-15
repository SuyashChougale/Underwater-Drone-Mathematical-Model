start_time = 0;
end_time = 10;
dt = 0.005;
times = start_time:dt:end_time;

N = numel(times);
x = [0; 0; 10];
xdot = zeros(3,1);
theta = zeros(3,1);
I = [IXX 0 0; 0 IYY 0; 0 0 IZZ];%Inertia matrix
%Some disturbance

deviation = 100;
thetadot = deg2rad(2*deviation*rand(3,1)-deviation);

for t = times
    i = input(t);
    omega = thetadot2omega(thetadot, theta);
    %a = acceleration(i, theta,xdot,m,g,k,kd);
    a = acceleration();
    %omegadot = angular_acceleration(i,omega, I,L,b,k);
    omegadot = angular_acceleration();

    omega = omega + dt * omegadot;
    thetadot = omega2thetadot(omega,theta);
    theta = theta + dt * thetadot;
    xdot = xdot + dt*a;
    x = x+ dt*xdot;

end

function omega=thetadot2omega()

end
function thetadot = omega2thetadot()

end
function a = acceleration()

end
function omegadot = angular_acceleration()

end
function t = thrust()

end
function t = torque()

end
