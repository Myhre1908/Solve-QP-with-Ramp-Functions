path = "C:\Users\jonas\OneDrive\Dokumenter\Matlab\5.klasse\Forprosjekt\QPramp\data\model1";
addpath("C:\Users\jonas\OneDrive\Dokumenter\Matlab\5.klasse\Forprosjekt\QPramp\mex")
input_folder = path;%"../data/model1";
initial_condition = 0;
timesteps = 20;

a = csvread(input_folder + "/a.csv"); %A_matrix
b = csvread(input_folder + "/b.csv"); %B-vector
x0 = csvread(input_folder + "/x0.csv"); %X0 values
invh = csvread(input_folder + "/invh.csv"); %A symmetric matrix
%invh = inv(C_A*Q_hat*C_A + R_hat)
w = csvread(input_folder + "/w.csv");
g = csvread(input_folder + "/g.csv");
s = csvread(input_folder + "/s.csv");
f = csvread(input_folder + "/f.csv");


% x0(1:2,:)
% s(1:2,:)
neg_g_invh_gt = -g*invh*g';
neg_s = -s;
neg_w = -w;
neg_g_invh = -g*invh;

nu = size(b, 2);
nx = size(a,2);
neg_invh_f = -invh*f;

% Transposing in the loop gives big slowdown
neg_g_invh_gt_t = neg_g_invh_gt';
neg_s_t = neg_s';
neg_w_t = neg_w';
neg_g_invh_t = neg_g_invh';
x0(1,:) = [-1,2.1];
x = x0(initial_condition+1,:)';

tic;
x_save = zeros(nx,timesteps);
x_save(:,1) = x;
u_save = zeros(nu,timesteps);
for i = 1:timesteps
    z = qpramp_solve(neg_g_invh_gt_t, neg_s_t, neg_w_t, neg_g_invh_t, x');
    u = z(1:nu)' + neg_invh_f(1:nu,:)*x;
    x = a*x + b*u;
    u_save(:,i) = u;
    x_save(:,i+1) = x;
end
t = toc;
fprintf("Total time for initial condition %d running %d timesteps: %.0f us\n", initial_condition, timesteps, t*1e6);
plot(x_save(1,:))
hold on 
plot(x_save(2,:))
clear mex; % Free memory used by solver
% plot(u_save)

