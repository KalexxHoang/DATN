%% Simulation RL for Mecanum Wheeled Mobile Robot
clc;clear;close all;
%% Time interval and simulation time
Step = 0.001;T_end = 150;
t = 0:Step:T_end;
%% Variables
q_r = cell(1,size(t,2));
q_r_dot = cell(1,size(t,2));

q = cell(1,size(t,2));
%q_dot = cell(1,size(t,2));

upsilon = cell(1,size(t,2)); % upsilon = v_q
upsilon_dot = cell(1,size(t,2));

upsilon_r = cell(1,size(t,2));
upsilon_r_dot = cell(1,size(t,2));

x = cell(1,size(t,2));
x_r = cell(1,size(t,2));
x_r_dot = cell(1,size(t,2));

e = cell(1,size(t,2));

z = cell(1,size(t,2));

W_c = cell(1,size(t,2));
W_a = cell(1,size(t,2));
gamma_upper = cell(1,size(t,2));

V_z = cell(1,size(t,2));
grad_V_z = cell(1,size(t,2));
tau = cell(1,size(t,2));
%% Parameter
m = 10; % Mass (kg)
l = 0.5; % Length of MWMR (m)
r = 0.05; % Wheeled radius (m)
I_theta = 5; I_phi = 0.1; % Rotational inertia (kg.m^2)
muy = 1; % Friction coefficients (N.m/rad/s)

ro = 0.1; % External disturbance
% ro = -0.1
Q = ones(6); % Positive-definite matrix
R = eye(3); % Positive-definite symmetric constant matrix
n_NN = 51; % Number of NN nodes
% The control gain
eta_c = 6.37;
eta_a1 = 7.8;
eta_a2 = 0.001;
eta_a = 0.001;
beta = 0.001;
gamma = 0.0001;
v = 0.002;
%% Initial Conditions
tau{1} = [13;14;-15];
q{1} = [1.8; 1.6; 1.4];
upsilon{1} = [0.3; 0.1; 0.2];
W_c{1} = 200*rand(n_NN,1);
W_a{1} = 200*rand(n_NN,1);
gamma_upper{1} = 10*eye(51);
%% System simulation
for i = 1:size(t,2)

if i<= 100000
        noise = sin(t(i))^2*cos(t(i)) + sin(2*t(i))^2*cos(0.1*t(i)) + sin(-1.2*t(i))^2*cos(0.5*t(i)) + sin(t(i))^5 + sin(1.2*t(i))^2 + cos(2.4*t(i))*sin(2.4*t(i))^3;
    else
        noise = 0;
end

m1 = ((4*m*r^2)/9) + (I_theta*r^2)/(9*l^2) + I_phi;
m2 = (I_theta*r^2)/(9*l^2) - (2*m*r^2)/9;
%M
M = [m1 m2 m2; m2 m1 m2; m2 m2 m1];
q_x = q{i}(1); q_y = q{i}(2); q_theta = q{i}(3);
% H^-1, H_cham
H_inv = [sin(q_theta + pi/3) -cos(q_theta + pi/3) -l; -sin(q_theta) cos(q_theta) -l; sin(q_theta-pi/3) -cos(q_theta-pi/3) -l];
dH_inv = [cos(q_theta + pi/3) sin(q_theta + pi/3) 0; -cos(q_theta) -sin(q_theta) 0; cos(q_theta - pi/3) sin(q_theta - pi/3) 0];
dH = (inv(H_inv))*dH_inv*(inv(H_inv));
%f_q
f_q = dH*H_inv - muy*(inv(H_inv))*(inv(M))*H_inv;
%g_q
g_q = r*(inv(H_inv))*(inv(M));
%fd(q,vq)

f_d_q_vq = -ro*(dH*H_inv*upsilon{i}+ g_q*tau{i});

%d(q,v_q)
d_q_vq = (inv(g_q))*f_d_q_vq;
% model
x{i} = [q{i}' upsilon{i}']'; % upsilon = v_q
f_x = [upsilon{i}' (f_q*upsilon{i})']';
g_x = [zeros(3) (g_q)']';
% G_z
G_z = [g_x; zeros(6,3)];
% g_qr
H_inv_qr = [sin(cos(0.5*t(i)) + pi/3) -cos(cos(0.5*t(i)) + pi/3) -l; -sin(cos(0.5*t(i))) cos(cos(0.5*t(i))) -l; sin(cos(0.5*t(i))-pi/3) -cos(cos(0.5*t(i))-pi/3) -l];
g_q_r = r*(inv(H_inv_qr))*(inv(M));
g_x_r = [zeros(3,3) g_q_r']';
g_cong = (inv(g_x_r'*g_x_r))*g_x_r';
%% reference trajectory
q_r{i} = [0.5*cos(2*t(i)) cos(t(i)) cos(0.5*t(i))]';
q_r_dot{i} = [-0.5*2*sin(2*t(i)) -sin(t(i)) -0.5*sin(0.5*t(i))]';
x_r{i} = [q_r{i}' q_r_dot{i}']';
% hr(xr) = đạo hàm x_r
x_r_dot{i} = [-sin(2*t(i)) -sin(t(i)) -0.5*sin(0.5*t(i)) -2*cos(2*t(i)) -cos(t(i)) -0.25*cos(0.5*t(i))]';
% f(xr)
dH_inv_qr = [cos(cos(0.5*t(i)) + pi/3) sin(cos(0.5*t(i)) + pi/3) 0; -cos(cos(0.5*t(i))) -sin(cos(0.5*t(i))) 0; cos(cos(0.5*t(i)) - pi/3) sin(cos(0.5*t(i)) - pi/3) 0];
dH_qr = (inv(H_inv_qr))*dH_inv_qr*(inv(H_inv_qr));
f_q_r = dH_qr*H_inv_qr - muy*(inv(H_inv_qr))*(inv(M))*H_inv_qr;
f_x_r = [q_r_dot{i}' (f_q_r*q_r_dot{i})']';

%% Caculate F_z
% tau_r
tau_r = g_cong*(x_r_dot{i}-f_x_r);
% F_z
F_z = [f_x-x_r_dot{i}+g_x*tau_r; x_r_dot{i}];
%% Tracking error
e{i} = x{i} - x_r{i};
z{i} = [e{i}' x_r{i}']';
%% Activation Function
% phi_z = zeros(51,1);
phi_z = (1/2)*[z{i}(1)^2 z{i}(2)^2 z{i}(3)^2 ... 
    z{i}(1)*z{i}(4) z{i}(1)*z{i}(5) z{i}(1)*z{i}(6) z{i}(2)*z{i}(4) z{i}(2)*z{i}(5) z{i}(2)*z{i}(6) ...  
    z{i}(3)*z{i}(4) z{i}(3)*z{i}(5) z{i}(3)*z{i}(6) (z{i}(1)^2)*(z{i}(2)^2) (z{i}(1)^2)*(z{i}(3)^2) ...  
    (z{i}(2)^2)*(z{i}(3)^2) (z{i}(1)^2)*(z{i}(7)^2) (z{i}(1)^2)*(z{i}(8)^2) (z{i}(1)^2)*(z{i}(9)^2) ... 
    (z{i}(1)^2)*(z{i}(10)^2) (z{i}(1)^2)*(z{i}(11)^2) (z{i}(1)^2)*(z{i}(12)^2) (z{i}(2)^2)*(z{i}(7)^2) ...  
    (z{i}(2)^2)*(z{i}(8)^2) (z{i}(2)^2)*(z{i}(9)^2) (z{i}(2)^2)*(z{i}(10)^2) (z{i}(2)^2)*(z{i}(11)^2) ... 
    (z{i}(2)^2)*(z{i}(12)^2) (z{i}(3)^2)*(z{i}(7)^2) (z{i}(3)^2)*(z{i}(8)^2) (z{i}(3)^2)*(z{i}(9)^2) ...  
    (z{i}(3)^2)*(z{i}(10)^2) (z{i}(3)^2)*(z{i}(11)^2) (z{i}(3)^2)*(z{i}(12)^2) (z{i}(4)^2)*(z{i}(7)^2) ... 
    (z{i}(4)^2)*(z{i}(8)^2) (z{i}(4)^2)*(z{i}(9)^2) (z{i}(4)^2)*(z{i}(10)^2) (z{i}(4)^2)*(z{i}(11)^2) ... 
    (z{i}(4)^2)*(z{i}(12)^2) (z{i}(5)^2)*(z{i}(7)^2) (z{i}(5)^2)*(z{i}(8)^2) (z{i}(5)^2)*(z{i}(9)^2) ... 
    (z{i}(5)^2)*(z{i}(10)^2) (z{i}(5)^2)*(z{i}(11)^2) (z{i}(5)^2)*(z{i}(12)^2) (z{i}(6)^2)*(z{i}(7)^2) ...
    (z{i}(6)^2)*(z{i}(8)^2) (z{i}(6)^2)*(z{i}(9)^2) (z{i}(6)^2)*(z{i}(10)^2) (z{i}(6)^2)*(z{i}(11)^2) ...
    (z{i}(6)^2)*(z{i}(12)^2)]';
% grad_phi_z = zeros(51,12);
grad_phi_z = (1/2)*[2*z{i}(1) 0 0 0 0 0 0 0 0 0 0 0; ...
    0 2*z{i}(2) 0 0 0 0 0 0 0 0 0 0; ...
    0 0 2*z{i}(3) 0 0 0 0 0 0 0 0 0; ...
    z{i}(4) 0 0 z{i}(1) 0 0 0 0 0 0 0 0; ...
    z{i}(5) 0 0 0 z{i}(1) 0 0 0 0 0 0 0; ...
    z{i}(6) 0 0 0 0 z{i}(1) 0 0 0 0 0 0; ...
    0 z{i}(4) 0 z{i}(2) 0 0 0 0 0 0 0 0; ...
    0 z{i}(5) 0 0 z{i}(2) 0 0 0 0 0 0 0; ...
    0 z{i}(6) 0 0 0 z{i}(2) 0 0 0 0 0 0; ...
    0 0 z{i}(4) z{i}(3) 0 0 0 0 0 0 0 0; ...
    0 0 z{i}(5) 0 z{i}(3) 0 0 0 0 0 0 0; ...
    0 0 z{i}(6) 0 0 z{i}(3) 0 0 0 0 0 0; ...
    2*z{i}(1)*(z{i}(2)^2) 2*(z{i}(1)^2)*z{i}(2) 0 0 0 0 0 0 0 0 0 0; ...
    2*z{i}(1)*(z{i}(3)^2) 0 2*(z{i}(1)^2)*z{i}(3) 0 0 0 0 0 0 0 0 0; ...
    0 2*z{i}(2)*(z{i}(3)^2) 2*(z{i}(2)^2)*z{i}(3) 0 0 0 0 0 0 0 0 0; ...
    2*z{i}(1)*(z{i}(7)^2) 0 0 0 0 0 2*(z{i}(1)^2)*z{i}(7) 0 0 0 0 0; ...
    2*z{i}(1)*(z{i}(8)^2) 0 0 0 0 0 0 2*(z{i}(1)^2)*z{i}(8) 0 0 0 0; ...
    2*z{i}(1)*(z{i}(9)^2) 0 0 0 0 0 0 0 2*(z{i}(1)^2)*z{i}(9) 0 0 0; ...
    2*z{i}(1)*(z{i}(10)^2) 0 0 0 0 0 0 0 0 2*(z{i}(1)^2)*z{i}(10) 0 0; ...
    2*z{i}(1)*(z{i}(11)^2) 0 0 0 0 0 0 0 0 0 2*(z{i}(1)^2)*z{i}(11) 0; ...
    2*z{i}(1)*(z{i}(12)^2) 0 0 0 0 0 0 0 0 0 0 2*(z{i}(1)^2)*z{i}(12); ...
    0 2*z{i}(2)*(z{i}(7)^2) 0 0 0 0 2*(z{i}(2)^2)*z{i}(7) 0 0 0 0 0; ...
    0 2*z{i}(2)*(z{i}(8)^2) 0 0 0 0 0 2*(z{i}(2)^2)*z{i}(8) 0 0 0 0; ...
    0 2*z{i}(2)*(z{i}(9)^2) 0 0 0 0 0 0 2*(z{i}(2)^2)*z{i}(9) 0 0 0; ...
    0 2*z{i}(2)*(z{i}(10)^2) 0 0 0 0 0 0 0 2*(z{i}(2)^2)*z{i}(10) 0 0; ...
    0 2*z{i}(2)*(z{i}(11)^2) 0 0 0 0 0 0 0 0 2*(z{i}(2)^2)*z{i}(11) 0; ...
    0 2*z{i}(2)*(z{i}(12)^2) 0 0 0 0 0 0 0 0 0 2*(z{i}(2)^2)*z{i}(12); ...
    0 0 2*z{i}(3)*(z{i}(7)^2) 0 0 0 2*(z{i}(3)^2)*z{i}(7) 0 0 0 0 0; ...
    0 0 2*z{i}(3)*(z{i}(8)^2) 0 0 0 0 2*(z{i}(3)^2)*z{i}(8) 0 0 0 0; ...
    0 0 2*z{i}(3)*(z{i}(9)^2) 0 0 0 0 0 2*(z{i}(3)^2)*z{i}(9) 0 0 0; ...
    0 0 2*z{i}(3)*(z{i}(10)^2) 0 0 0 0 0 0 2*(z{i}(3)^2)*z{i}(10) 0 0; ...
    0 0 2*z{i}(3)*(z{i}(11)^2) 0 0 0 0 0 0 0 2*(z{i}(3)^2)*z{i}(11) 0; ...
    0 0 2*z{i}(3)*(z{i}(12)^2) 0 0 0 0 0 0 0 0 2*(z{i}(3)^2)*z{i}(12); ...
    0 0 0 2*z{i}(4)*(z{i}(7)^2) 0 0 2*(z{i}(4)^2)*z{i}(7) 0 0 0 0 0; ...
    0 0 0 2*z{i}(4)*(z{i}(8)^2) 0 0 0 2*(z{i}(4)^2)*z{i}(8) 0 0 0 0; ...
    0 0 0 2*z{i}(4)*(z{i}(9)^2) 0 0 0 0 2*(z{i}(4)^2)*z{i}(9) 0 0 0; ...
    0 0 0 2*z{i}(4)*(z{i}(10)^2) 0 0 0 0 0 2*(z{i}(4)^2)*z{i}(10) 0 0; ...
    0 0 0 2*z{i}(4)*(z{i}(11)^2) 0 0 0 0 0 0 2*(z{i}(4)^2)*z{i}(11) 0; ...
    0 0 0 2*z{i}(4)*(z{i}(12)^2) 0 0 0 0 0 0 0 2*(z{i}(4)^2)*z{i}(12); ...
    0 0 0 0 2*z{i}(5)*(z{i}(7)^2) 0 2*(z{i}(5)^2)*z{i}(7) 0 0 0 0 0; ...
    0 0 0 0 2*z{i}(5)*(z{i}(8)^2) 0 0 2*(z{i}(5)^2)*z{i}(8) 0 0 0 0; ...
    0 0 0 0 2*z{i}(5)*(z{i}(9)^2) 0 0 0 2*(z{i}(5)^2)*z{i}(9) 0 0 0; ...
    0 0 0 0 2*z{i}(5)*(z{i}(10)^2) 0 0 0 0 2*(z{i}(5)^2)*z{i}(10) 0 0; ...
    0 0 0 0 2*z{i}(5)*(z{i}(11)^2) 0 0 0 0 0 2*(z{i}(5)^2)*z{i}(11) 0; ...
    0 0 0 0 2*z{i}(5)*(z{i}(12)^2) 0 0 0 0 0 0 2*(z{i}(5)^2)*z{i}(12); ...
    0 0 0 0 0 2*z{i}(6)*(z{i}(7)^2) 2*(z{i}(6)^2)*z{i}(7) 0 0 0 0 0; ...
    0 0 0 0 0 2*z{i}(6)*(z{i}(8)^2) 0 2*(z{i}(6)^2)*z{i}(8) 0 0 0 0; ...
    0 0 0 0 0 2*z{i}(6)*(z{i}(9)^2) 0 0 2*(z{i}(6)^2)*z{i}(9) 0 0 0; ...
    0 0 0 0 0 2*z{i}(6)*(z{i}(10)^2) 0 0 0 2*(z{i}(6)^2)*z{i}(10) 0 0; ...
    0 0 0 0 0 2*z{i}(6)*(z{i}(11)^2) 0 0 0 0 2*(z{i}(6)^2)*z{i}(11) 0; ...
    0 0 0 0 0 2*z{i}(6)*(z{i}(12)^2) 0 0 0 0 0 2*(z{i}(6)^2)*z{i}(12)];
%% Estimated controlling performance
V_z{i} = W_c{i}'*phi_z;
grad_V_z{i} = grad_phi_z'*W_c{i}; % Xem lại đạo hàm
%% Controller
u = -1/2*(inv(R))*G_z'*grad_phi_z'*W_a{i} + Step*noise;
%u = -1/2*(inv(R))*G_z'*grad_phi_z'*W_a{i};
tau{i+1} = u + g_cong*(x_r_dot{i} - f_x_r);


%% Update stage
F_q = upsilon{i};
F_upsilon = f_q*upsilon{i} + g_q*tau{i} + g_q*d_q_vq;

if i==size(t,2)
   break;
end
upsilon{i+1} = upsilon{i} + Step*F_upsilon;
q{i+1} = q{i} + Step*F_q;
% upsilon{i+1} = upsilon{i} + Step*F_upsilon + Step*noise;
% q{i+1} = q{i} + Step*F_q + Step*noise;
%% Update NN Weight
% W_c
sigma = grad_phi_z*(F_z+G_z*u)-gamma*phi_z;
D_1 = grad_phi_z*G_z*(inv(R))*G_z'*grad_phi_z';
delta = (0.01*(x{i}(1)^2 + x{i}(2)^2 +x{i}(3)^2 + x{i}(4)^2 + x{i}(5)^2 + x{i}(6)^2)*norm(g_x))^2 + z{i}'*[Q zeros(6); zeros(6) zeros(6)]*z{i} + 1/4*W_a{i}'*D_1*W_a{i} + W_c{i}'*(sigma);
%delta = z{i}'*[Q zeros(6); zeros(6) zeros(6)]*z{i} + 1/4*W_a{i}'*D_1*W_a{i} + W_c{i}'*(sigma);
W_c{i+1} = W_c{i} + Step*(-eta_c*gamma_upper{i}*sigma/(1+v*sigma'*gamma_upper{i}*sigma)*delta);
gamma_upper{i+1} = gamma_upper{i} - eta_c*(-beta*gamma_upper{i}+gamma_upper{i}*(sigma*sigma')/(1+v*sigma'*gamma_upper{i}*sigma)*gamma_upper{i});

% W_a
sigma_ngang = sigma/(1+v*sigma'*gamma_upper{i}*sigma)^(1/2);
W_a{i+1} = W_a{i} + Step*(-eta_a1*(W_a{i} - W_c{i})-eta_a2*W_a{i}+eta_a/4*D_1*W_a{i}*sigma_ngang'/((1+v*sigma'*gamma_upper{i}*sigma)^(1/2))*W_c{i});
end

%% Plot

q = cell2mat(q);
q_r = cell2mat(q_r);

figure;
subplot(3,1,1);
plot(t,q(1,:));
hold on
plot(t,q_r(1,:));
title('X position');
hold off
subplot(3,1,2);
plot(t,q(2,:));
hold on
plot(t,q_r(2,:));
title('Y position');
hold off
subplot(3,1,3);
plot(t,q(3,:));
hold on
plot(t,q_r(3,:));
title('Heading angle');
hold off