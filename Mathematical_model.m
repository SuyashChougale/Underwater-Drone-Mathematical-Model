%time initiation
start_time = 0;
end_time = 10;
dt = 0.005;
times = start_time:dt:end_time;

N = numel(times);%number of instances
x = [0; 0; 10];
xdot = zeros(3,1);
theta = zeros(3,1);

% Initialize storage for plotting
x_store = zeros(3, N);
xdot_store = zeros(3, N);
theta_store = zeros(3, N);
thetadot_store = zeros(3, N);
a_store = zeros(3,N);
% Some disturbance
deviation = 0;
thetadot = deg2rad(2*deviation*rand(3,1)-deviation);


m = 1.04; % Mass of the quadcopter (kg)
g = 9.81; % Acceleration due to gravity (m/s^2)
I=[0.01122 0 0; 0 0.01122 0;0 0 0.01122];

for idx = 1:N
    S = 0.006;
    if x(3)>0
        Cd = 0.63;
        rho_medium = 1.135;
    else
        % The drag coefficient is found to be about 1.8 at Reynolds number of 1000, which is at a flow velocity of 0.001m/s.
        Cd = 1.8;
        rho_medium = 997;
    end
    kd = Cd*rho_medium*S/2;
    L = 0.23;
    R = 0.127;
    if x(3)>0
        Ct = 0.137;
        rho_medium = 1.135;
    else
        % The drag coefficient is found to be about 1.8 at Reynolds number of 1000, which is at a flow velocity of 0.001m/s.
        Ct = 1.8;% look value of thrust coefficinets in book
        rho_medium = 997;
    end
    k = Ct*rho_medium*R*R*R*R/2/2/pi/pi;
    
    if x(3)>0
        Ch = 0.0092;
        rho_medium = 1.135;
    else
        % The drag coefficient is found to be about 1.8 at Reynolds number of 1000, which is at a flow velocity of 0.001m/s.
        Ch = 0.0092;% look value of thrust coefficinets in book
        rho_medium = 997;
    end
    b = Ch*rho_medium*R*R*R*R*R/2/2/pi/pi;
    t = times(idx);
    i = input(t); % Assuming input(t) is a function providing the inputs at time t
    omega = thetadot2omega(thetadot, theta);
    a = acceleration(i, theta, xdot, m, g, k, kd,x(3));
    
    omegadot = angular_acceleration(i, omega, I, L, b, k);
    
    omega = omega + dt * omegadot;
    thetadot = omega2thetadot(omega, theta);
    theta = theta + dt * thetadot;
    xdot = xdot + dt * a;
    x = x + dt * xdot;
    
    % Store values
    x_store(:, idx) = x;
    xdot_store(:, idx) = xdot;
    theta_store(:, idx) = theta;
    thetadot_store(:, idx) = thetadot;
    a_store(:, idx) = a;
end

% Plot results
figure;
subplot(5,1,1);
plot(times, x_store);
title('Position (x) over Time');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');

subplot(5,1,2);
plot(times, xdot_store);
title('Velocity (xdot) over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('xdot_x', 'xdot_y', 'xdot_z');

subplot(5,1,4);
plot(times, theta_store);
title('Orientation (theta) over Time');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('theta_x', 'theta_y', 'theta_z');

subplot(5,1,5);
plot(times, thetadot_store);
title('Angular Velocity (thetadot) over Time');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('thetadot_x', 'thetadot_y', 'thetadot_z');

subplot(5,1,3);
plot(times, a_store);
title('Acceleration over Time');
xlabel('Time (s)');
ylabel('Acceleration (m^2/s)');
legend('a_x', 'a_y', 'a_z');

% Functions (unchanged from previous version)
function T = thrust(inputs, k)
    T = [0; 0; k * sum(inputs)];
end

function tau = torque(inputs, L, b, k)
    tau = [L * k * (inputs(1) - inputs(3)); 
           L * k * (inputs(2) - inputs(4)); 
           b * (inputs(1) - inputs(2) + inputs(3) - inputs(4))];
end

function omega = thetadot2omega(thetadot, theta)
    phi = theta(1);
    theta_angle = theta(2);
    psi = theta(3);
    
    R = [1, 0, -sin(theta_angle);
         0, cos(phi), cos(theta_angle) * sin(phi);
         0, -sin(phi), cos(theta_angle) * cos(phi)];
     
    omega = R * thetadot;
end

function thetadot = omega2thetadot(omega, theta)
    phi = theta(1);
    theta_angle = theta(2);
    
    R_inv = [1, sin(phi) * tan(theta_angle), cos(phi) * tan(theta_angle);
             0, cos(phi), -sin(phi);
             0, sin(phi) / cos(theta_angle), cos(phi) / cos(theta_angle)];
         
    thetadot = R_inv * omega;
end

function a = acceleration(inputs, angles, xdot, m, g, k, kd,x)
    rho = 997;
    V =  6.9656e-04;
    if x>0
      gravity = [0; 0; -g];
    else
      gravity = [0; 0; rho*V*g/m-g];
    end
    % gravity = [0; 0; -g];
    R = rotation(angles);
    T = R * thrust(inputs, k);
    Fd = -kd * xdot;
    a = gravity + (1/m) * T + Fd;
end

function R = rotation(angles)
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);
    
    R = [cos(psi) * cos(theta), cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi), cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
         sin(psi) * cos(theta), sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi), sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
         -sin(theta), cos(theta) * sin(phi), cos(theta) * cos(phi)];
end

function omegadot = angular_acceleration(inputs, omega, I, L, b, k)
    tau = torque(inputs, L, b, k);
    omegadot = inv(I) * (tau - cross(omega, I * omega));
end

function i = input(t)
    % Placeholder for the input function. Replace with actual input logic.
    a = 71;
    b = 71;
    c = 71;
    d = 71;
    i = [a; b; c; d]; % Example input, same for all motors

end
