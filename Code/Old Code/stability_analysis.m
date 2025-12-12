
d = 4e-3;
ke = 5;

de = 0.1; %this combination induces a hopf bifurcation
va = 5; %was 4.8e-7
da = 0.5;
Nm = 1.2; % might be e-9
k = 0.9;

N = de*da*k/(d*ke*va);
E = k*de/(ke*d);
A = k/d;
J = [k*(1 - 2*N/Nm)-d*E, -d*N, 0;
     0, -de, ke;
     va, 0, -da];
% 
C = zeros(3,3);
C(1,1) = 1;
D = zeros(3,3);
% system = ss(J, -0.5*eye(3), C, D);
% bode(system)
% 
% eig(J)

% B = -0.5*eye(3);

% eq1 = (k*(1-2*N/Nm) - d*E - L)*((de + L)*(da+L)) - va*ke == 0;
% eq2 = (k-d*E)*N == 0;
% eq3 = ke*A - de*E == 0;
% eq4 = va*N - da*A == 0;

% syms k d de da va ke N E A L Nm
% eq2 = A == k*de / (ke*d)
% eq3 = E == k/d
% eq4 = N == de*da*k / (d*ke*va)
% % 

% solve(eq1, L)


%% Plotting Induced Identification:
clear
clc 
close all

d = 100e-3;
ke = 5;

de = 0.1; %this combination induces a hopf bifurcation
va = 5; %was 4.8e-7
da = 0.5;
Nm = 1.2; % might be e-9
k = 0.9;

N = de*da*k/(d*ke*va);
E = k*de/(ke*d);
A = k/d;
J = [k*(1 - 2*N/Nm)-d*E, -d*N, 0;
     0, -de, ke;
     va, 0, -da];
eig(J)
init = [1.4 1 3];
tspan = [0 1000];

[t, x] = ode45(@dxdt, tspan, init);
nexttile;
plot(t,x(:,1))
nexttile;
plot(t,x(:,2))
nexttile
plot(t,x(:,3))

function derivs=dxdt(t,x)

d = 100e-3;
ke = 5;

de = .1; %this combination induces a hopf bifurcation
va = 5; %was 4.8e-7
da = 0.5;
Nm = 1.2; % might be e-9
k = 0.9;

derivs = zeros(3,1);
derivs(1) = k*x(1)*(1 - x(1)/Nm)-d*x(2)*x(1);
derivs(2) = ke*x(3) - de*x(2);
derivs(3) = va*x(1) - da*x(3);

end