%FYE take home 2019 retake Q2 part (ii)
% solve 1D Heat eqn with FTCS method.
% then compare numerical results vs analytic soluntions
% u_t = u_xx

% Jianhong Chen
% 09-20-2019

clear all
clc

% problem setup parameters
dx = 1/32; % given
dt = (dx^2)/2; % compute the stable time step 
t0 = 0;
T = 5; % Terminal time that is long enough to reach steady state
a = 0;
b = 1;
Nx = (b-a)/dx - 1; %internal grid
Nt = (T-t0)/dt; % time step point
r = dt/(dx^2); 

% IC functions
f1 = @(x) 2*x;
f2 = @(x) 2*(1-x);

% analytic soln funciton
F = @(n,x,t) 8/(n*pi).^2 .*sin(n*pi/2).*sin(n*pi*x).*exp(-(n*pi).^2.*t); 

x = a + (1:Nx)*dx;%construct x-step vector
x = [a;x.';b]; % transpose of A in matlab A.'

% allocate matrix-vector arrays
u = zeros(Nx+2, Nt);

% BC
u(1, :) = 0;
u(Nx+2, :) = 0;
% IC
u(2:17, 1) =f1(x(2:17));
u(17:32,1) = f2(x(17:32));

% important!!! the whole Heat euqantion in space has t obe solve for each
% time step. 
for n = 1:Nt % time has to be in the outer loop
    for i = 2:Nx+1 % solve all space first for every timeframe
        u(i,n+1) = u(i,n) + r*(u(i+1, n) - 2*u(i,n) + u(i-1, n));
    end
end

% find t_max by computing the temporal change of the numerical soln
n_ss = 1; % steady state step 
sigma = 10^-6;
u_res = 1; %tempolar change residue

while u_res > sigma
    u_res = sum(abs(u(:,n_ss+1) - u(:,n_ss)))/(dt*(Nx+1)); 
    n_ss = n_ss + 1;
end

t_max = n_ss*dt; %  compute time that take to reach steady state

T_save = [0.25, 0.5, 0.75, 1]*t_max;
n_save = round(T_save/dt);

% now compute analytic soln as the exact soln
u_ana = zeros(Nx+2, length(T_save));
n_ana = 1:30; % only sum up the first 30 terms

for t = 1:length(T_save) % time-loop
    u_dummy = zeros(length(n_ana),1);
    for i = 1:length(x) % x-loop
        for n = n_ana % compute each term for n
            u_dummy(n) = F(n,x(i),T_save(t));            
        end
        u_ana(i,t) = sum(u_dummy);% sum up all n-terms
    end
end

% save workspace variables
save dat_heat_FTCS

% 
figure(1)
clf
hold on
plot(x, u(:, n_save(1)), 'ro', 'linewidth',2)
plot(x, u(:, n_save(2)), 'bs','linewidth',2)
plot(x, u(:, n_save(3)), 'g+', 'linewidth',2)
plot(x, u(:, n_save(4)), 'cx', 'linewidth',2)

plot(x, u_ana(:,1), 'r-', 'linewidth',2)
plot(x, u_ana(:,2), 'b-', 'linewidth',2)
plot(x, u_ana(:,3), 'g-', 'linewidth',2)
plot(x, u_ana(:,4), 'c-', 'linewidth',2)

legend('num-t/tmax=0.25', 'num-t/tmax=0.5', ...
    'num-t/tmax=0.75', 'num-t/tmax=1', ...
    "ana-t/tmax = 0.25", "ana-t/tmax = 0.5", ...
    "ana-t/tmax = 0.75", "ana-t/tmax = 1")
xlabel("x")
ylabel("u(x,t)")
title("1D Heat with FTCS at Different Times $r = \frac{1}{2}$", ... 
    'interpreter', 'latex')




