%FYE take home 2019 retake Q2 part (ii)
% solve 1D Heat eqn with BTCS method.
% then compare numerical results vs analytic soluntions
% u_t = u_xx

% Jianhong Chen
% 09-20-2019

clear all
clc

% load the same parameters from previous solution
load dat_heat_FTCS

% allocate matrix-vector arrays
u_BT = zeros(Nx+2, n_save(end));

% BC
u_BT(1, :) = 0;
u_BT(Nx+2, :) = 0;
% IC
u_BT(2:17, 1) =f1(x(2:17));
u_BT(17:32,1) = f2(x(17:32));

% BTCS method has the form 
% u(n+1) = u(n) + r*(Ax*u(n+1)
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
% rearrange BTCS into Ax = b form
% (I-r*Ax)*u(n+1) = u(n)
% thus A = (I - r*Ax)
A = speye(Nx,Nx) - r*Ax;
clear Ax

for n = 1:n_save(end) % time-loop
    % only update internal points 
    % in this case it doesn't matter, becasue BC are zeros
    b = u_BT(2:Nx+1,n);
    %u_BT(2:Nx+1,n+1) = A\b;
    u_BT(2:Nx+1,n+1) = GaussianElimination(A, b); % solve Ax = b by GE
end
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

plot(x, u_BT(:,n_save(1)), 'r--', 'linewidth',2)
plot(x, u_BT(:,n_save(2)), 'b--', 'linewidth',2)
plot(x, u_BT(:,n_save(3)), 'g--', 'linewidth',2)
plot(x, u_BT(:,n_save(4)), 'c--', 'linewidth',2)

legend('FT-t/tmax=0.25', 'FT-t/tmax=0.5', ...
    'FT-t/tmax=0.75', 'FT-t/tmax=1', ...
    "ana-t/tmax = 0.25", "ana-t/tmax = 0.5", ...
    "ana-t/tmax = 0.75", "ana-t/tmax = 1", ...
    'BT-t/tmax=0.25', 'BT-t/tmax=0.5', ...
    'BT-t/tmax=0.75', 'BT-t/tmax=1')
xlabel("x")
ylabel("u(x,t)")
title("1D Heat at Different Times $r = \frac{1}{2}$", ... 
    'interpreter', 'latex')





