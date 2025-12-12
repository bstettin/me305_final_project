clear
clc
close all

d = 4e-3;
ke = 5;
de = 2;
va = 4.8e-7;
da = 0.5;
Nm = 1.2; 
k = 0.9;
Ae = k*de / (d*ke);
% 
% Ee = k / d;
% Ae = de/ke*Ee;
% Ne = da*de*k/(va*ke*d);

x0 = [pi, 0.25];
x = fsolve(@hill_zeros, x0);
disp(['Omega = ' num2str(x(1)) ' rad'])
disp(['  Tau = ' num2str(x(2)) ' hrs'])






%%
% 2D system:
delays = [0.2];
history = [0.1 0];
tspan = [0, 525];
sol = dde23(@derivs2D, delays, history, tspan);

plot2D(sol,history, x(1),x(2))

function roots = hill_zeros(x)
    %initialize model parameters:
    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = 0.5;
    Nm = 0.5; 
    k = 2;
    Ae = k*de / (d*ke);
    kh = 0.1;


    omega = x(1);
    tau = x(2);
    roots = zeros(2,1);

    roots(1) = da*k - da/kh*cos(omega*tau) - omega/kh*sin(omega*tau)+omega^2;
    roots(2) = da/kh*sin(omega*tau) - omega*da + omega*k - omega/kh*cos(omega*tau);
end

function dydt = derivs2D(t, y, Z)
    N_tau = Z(1,1);

    d = 4e-3;
    ke = 5;

    de = 2;
    va = 4.8e-7;
    da = 0.5;
    Nm = 0.5; % might be e-9
    k = 2; %0.8
    kh = 0.1;

    % Ae = k*de / (d*ke);

    % y = [N E A]
    N = y(1);
    A = y(2);
    E = ke/de*A;

    dNdt = k * N - k*N^2 / Nm - d*E*N - N_tau / (kh + N_tau);
    dAdt = va*N- da*A;

    dydt = [dNdt;
            dAdt];
end