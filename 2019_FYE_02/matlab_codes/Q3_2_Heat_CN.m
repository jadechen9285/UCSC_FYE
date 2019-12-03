%FYE take home 2019 retake Q3 part (ii)
% solve 1D Heat eqn with Crank-Nicolson method.
% then compare numerical results vs analytic soluntions
% u_t = u_xx

% Jianhong Chen
% 09-20-2019

clear all
clc

% problem setup parameters
N = 32;
dx = 1/N; % given
dt = 5*10^-4; % compute the stable time step 
t0 = 0;
T = 5;
a = 0;
b = 1;
Nx = (b-a)/dx - 1; %internal grid
Nt = round((T-t0)/dt);
r = dt/(dx^2);

% IC functions
f1 = @(x) 2*x;
f2 = @(x) 2*(1-x);

% analytic soln funciton
F = @(n,x,t) 8/(n*pi).^2 .*sin(n*pi/2).*sin(n*pi*x).*exp(-(n*pi).^2.*t); 

x = a + (1:Nx)*dx;
x = [a;x.';b]; % transpose of A in matlab A.'

u_CN = zeros(Nx+2, Nt);
% BC
u_CN(1, :) = 0;
u_CN(Nx+2, :) = 0;
% IC
u_CN(2:17, 1) =f1(x(2:17));
u_CN(17:32,1) = f2(x(17:32));

% CTCS method has the form 
% u(n+1) = u(n) + r/2*(Ax(u(n)) +  r/2*(Ax*u(n+1))
% due to zero BC, b-vector is zero in this case

% contruct sparse matrix Ax
Ia=zeros(1,3*Nx);    % for storing row indices of non-zero entries
Ja=zeros(1,3*Nx);    % for storing column indices of non-zero entries
Sa=zeros(1,3*Nx);    % for storing values of non-zero entries
for i=1:Nx
  Ia(3*(i-1)+[1:3])=[i,i,i];
  Ja(3*(i-1)+[1:3])=[i-1,i,i+1];
  Sa(3*(i-1)+[1:3])=[1, -2, 1];
end
ind=find( (Ja>0)&(Ja<Nx+1) );
Ax=sparse(Ia(ind),Ja(ind),Sa(ind),Nx,Nx);
% rearrange CN method into Ax = b form
% (I-r/2*Ax)*u(n+1) = (I+r/2*Ax)*u(n)
% thus A = (I - r*Ax)
A = speye(Nx,Nx) - r/2*Ax;

% impementing Jacobi iteration scheme

%construct D and R matrix
D = diag(diag(A)); % only the diagonal elements of A
R = A - D; % the rest elements except the diagonal ones of A
% inverse of diagonal matrix is simply the reciprocal 
D_inv = diag(1./diag(A));
% compute J-matrix
J = -D_inv * R;

% iterate the Jacobi scheme forward in time
N_Jacobi = zeros(Nt,1); %record the step for convergent for Jacobi method
for n = 1:Nt
    b = (speye(Nx,Nx)+r/2*Ax)*u_CN(2:Nx+1, n);
    c = D_inv *b;
    sigma = 10^-6;
    X = zeros(Nx, 100);
    J_res = 1; % Jacobi iteration error
    k = 1;
    X(:,1) = 1; % initial guess of X for Jacobi scheme
    while J_res > sigma
        X(:,k+1) = J*X(:,k) + c;
        J_res = 1/(dt*N)*sum(abs(X(:,k+1) - X(:,k)));
        k = k +1;
    end
    N_Jacobi(n) = k; %save max Jacobi iteration step
    u_CN(2:Nx+1, n+1) = X(:,k);
    % check if steady state is reached
    u_res = 1/(dt*N)*sum(abs(u_CN(:,n+1) - u_CN(:,n))); 
    if u_res < sigma
        n_ss = n;
        break % stop the time-loop once reach steady state
    end
        
end

t_max = n_ss*dt;
T_save = [0.25, 0.5, 0.75, 1]*t_max;
n_save = round(T_save/dt);

% now compute analytic soln for as the benchmark comparison
u_ana = zeros(Nx+2, length(T_save));
n_ana = 1:30; % only sum up the first 30 terms

for t = 1:length(T_save)
    u_temp = zeros(length(n_ana),1);
    for i = 1:length(x)
        for n = n_ana
            u_temp(n) = F(n,x(i),T_save(t));            
        end
        u_ana(i,t) = sum(u_temp);
    end
end

% compute error of CN with Jacobi scheme
error = sum(abs(u_CN(:,n_save(end)) - u_ana(:,4)))/N;
disp("the error at t_max is ")
disp(error)


figure(1)
clf
hold on

plot(x, u_CN(:, n_save(1)), 'ro', 'linewidth',2)
plot(x, u_CN(:, n_save(2)), 'bs','linewidth',2)
plot(x, u_CN(:, n_save(3)), 'g+', 'linewidth',2)
plot(x, u_CN(:, n_save(4)), 'cx', 'linewidth',2)

plot(x, u_ana(:,1), 'r-', 'linewidth',2)
plot(x, u_ana(:,2), 'b-', 'linewidth',2)
plot(x, u_ana(:,3), 'g-', 'linewidth',2)
plot(x, u_ana(:,4), 'c-', 'linewidth',2)

legend('CN-t/tmax=0.25', 'CN-t/tmax=0.5', ...
    'CN-t/tmax=0.75', 'CN-t/tmax=1', ...
    "ana-t/tmax = 0.25", "ana-t/tmax = 0.5", ...
    "ana-t/tmax = 0.75", "ana-t/tmax = 1")

xlabel("x")
ylabel("u(x,t)")
title("1D Heat at Different Times with $\Delta t = 5* 10^{-4}$", ... 
    'interpreter', 'latex')










