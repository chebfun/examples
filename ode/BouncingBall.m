%% A Bouncing Ball
% Filomena Di Tommaso, February 26th, 2013

%%
% (Chebfun example ode/BouncingBall.m)
% [Tags: #linearODE, #Ball, #physics]

%%
% This example simulates a bouncing ball. The trajectory of the ball is a
% piecewise quadratic chebfun.

%%
% Initialise:
clc, clear all, close all
trajectory = chebfun;

%%
% Parameters:
g = 9.81;                   % Force of gravity
x0 = 0;                     % Initial point
y0 = 0;
v0 = 30;                    % Initial velocity
theta(1) = pi/4;            % Starting angle
k = 0.9;                    % Coefficient of rebound
r = 1;                      % Radius of the ball

i = 1;                      % Counter of the number of parabolas
v_0x = v0*cos(theta);       % Initial velocity in the x direction
v_0y(i) = v0*sin(theta);    % Initial velocity in the y direction

%% 
% The main loop:

while ( theta > 0.15 )
    
    % Time at which the ball impacts the ground:
    t_ground(i) = (v_0y(i)+sqrt(v_0y(i)^2+2*g*y0))/g;
    
    % x coordinate of the impact point:
    x_ground(i) = x0+v_0x*t_ground(i);
    
    % Motion law;
    x = chebfun('x', [x0, x_ground(i)]);
    y = y0 + (-g/(2*v_0x^2))*(x-x0).^2 + (v_0y(i)/v_0x)*(x-x0);
    
    % Piecewise chebfun:
    trajectory = [trajectory ; y];
    
    % Velocity at the impact point:
    v_fy=v_0y(i)-g*t_ground(i);
    
    % The y component of the velocity decrease according to the coefficient of
    % rebound k:
    v_0y(i+1)=-k*v_fy;
    
    x0 = x_ground(i);
    y0 = 0;
    
    % Corner rebound ball:
    theta(i) = atan(v_0y(i+1)/v_0x);
    
    i = i+1;
end

%%
% The final part of the trajectory is horizontal:
final_part = chebfun(@(y) 0*y);
trajectory = [trajectory ; final_part];

%%
% I shift the trajectory vertically by a quantity equal to the radius r in order
% to visualize the impact with the ground:
trajectory = trajectory + r;

%%
% Parameter to plot the circle:
n = 50;
t = linspace(0, 2*pi, n); 
t_ground = [0, t_ground, 5];
x_ground = [0, x_ground, 900];
% pause(2)

c = 1;
ks = 120;   % Coefficient of elasticity of a ball
m = 0.313;  % Weight of the ball in gr

%%
% Plot the result nicely:
axis([0, 900 , 0, 25])

% Since I can't set axis equal to see the ball I must increase the eccentricity.
a = 25;
b = 1;
M = [];
for j = 1:i
    % Temporal interval
    time = linspace(0, t_ground(j+1), n);
    xt = x_ground(j) + v_0x*time;
    ballx = xt(1) + a*r*cos(t);
    bally = trajectory(xt(1)) + b*r*sin(t);
    h = patch(ballx, bally, 'g'); hold on
    g = plot(ballx(1), bally(1), '*r'); % Point on the ball
    drawnow
%     M = [M getframe()];
    for k = 2:n
        if ( b < 1 )
            b = b + ((1-b)/(i)*j);
        else
            b = 1;
        end
        ballx = xt(k) + a*r*cos(t);
        bally = trajectory(xt(k)) + b*r*sin(t);
        set(h, 'XData', ballx, 'YData', bally)
        set(g, 'XData', ballx(n-k+1), 'YData', bally(n-k+1))
        drawnow
%         M = [M getframe()];
    end
    b = r - (1/2*m*v_0y(j)^2)/ks;
    ballx = xt(k) + a*r*cos(t);
    bally = trajectory(xt(k)) + b*r*sin(t);
    set(h, 'XData', ballx, 'YData', bally)
    set(g, 'XData', ballx(n-k+1), 'YData', bally(n-k+1))
    drawnow
%     M = [M getframe()];
    
    set(h, 'XData', ballx, 'YData', bally, 'Visible', 'off')
    set(g, 'XData', ballx(n-k+1), 'YData', bally(n-k+1), 'Visible', 'off')
    
    
end


%%
% We now write the images to an animated .gif!

% for j = 1:5:numel(M)
%     im = frame2im(M(j));
%     [imind,cm] = rgb2ind(im,16);
%     if j == 1;
%         imwrite(imind, cm, 'html/BouncingBall.gif', 'gif', ...
%             'Loopcount', inf, 'DelayTime', 1e-5);
%     else
%         imwrite(imind, cm, 'html/BouncingBall.gif', 'gif', ...
%             'WriteMode', 'append', 'DelayTime', 1e-5);
%     end
% end
