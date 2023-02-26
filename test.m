%% Simulation RL for Mecanum Wheeled Mobile Robot
    clc;clear;close all;
%% Time interval and simulation time
    Step = 0.001;T_end = 150;
    t = 0:Step:T_end;
     q_r = cell(1,size(t,2));
    q_r_dot = cell(1,size(t,2));
     x_r = cell(1,size(t,2));
    x_r_dot = cell(1,size(t,2));

for i = 1:size(t,2)
    q_r{i} = [0.5*cos(2*t(i)) -0.5*sin(2*t(i)) cos(0.5*t(i))]';
    q_r_dot{i} = [-sin(2*t(i)) -cos(2*t(i)) -0.5*sin(0.5*t(i))]';
    x_r{i} = [q_r{i}' q_r_dot{i}']';
    % dxr = hr(xr)
    x_r_dot{i} = [-sin(2*t(i)) -cos(2*t(i)) -0.5*sin(0.5*t(i)) -2*cos(2*t(i)) 2*sin(2*t(i)) -0.25*cos(0.5*t(i))]';
end
    q_r = cell2mat(q_r);
    plot(q_r(1,:),q_r(2,:));