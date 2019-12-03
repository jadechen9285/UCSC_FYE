%FYE take home 2019 retake Q2 part (ii)
% solve 1D Heat eqn with BTCS method.
% then analysis number of steps it takes to reach steady state vs vary 
% time step dt.
% u_t = u_xx

% Jianhong Chen
% 09-20-2019

clear all
clc

% load the same parameters from previous solution

dx = 1/32; % given
C = 10; %constant choice for time step
dt = (dx^2)*C/2; % compute the stable time step 
t0 = 0;
T = 5;
a = 0;
b = 1;
Nx = (b-a)/dx - 1; %internal grid
Nt = (T-t0)/dt;
r = dt/(dx^2);

% IC functions
f1 = @(x) 2*x;
f2 = @(x) 2*(1-x);

% analytic soln funciton
F = @(n,x,t) 8/(n*pi).^2 .*sin(n*pi/2).*sin(n*pi*x).*exp(-(n*pi).^2.*t); 

x = a + (1:Nx)*dx;
x = [a;x.';b]; % transpose of A in matlab A.'

% allocate matrix-vector arrays
u_BT = zeros(Nx+2, Nt);

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
A = speye(Nx,Nx) - r(1)*Ax;
clear Ax

sigma = 10^-6;

for n = 1:Nt % time-loop
    % only update internal points 
    % in this case it doesn't matter, becasue BC are zeros
    b = u_BT(2:Nx+1,n);
    %u_BT(2:Nx+1,n+1) = A\b;
    u_BT(2:Nx+1,n+1) = GaussianElimination(A, b); % solve Ax = b by GE
    u_res = 1/(dt*Nx+1)* sum(abs(u_BT(:,n+1) - u_BT(:,n)));
    if u_res < sigma
        n_ss = n;
        break % break the time loop once reach steady state
    end
end

t_max = n_ss*dt;







