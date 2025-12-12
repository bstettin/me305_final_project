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

Ee = k / d
Ae = de/ke*Ee
Ne = da*de*k/(va*ke*d)
% 
% Ae = k*de / (d*ke);
% Ee = ke*Ae / de;


% x0 = [pi/10, 10];
% x = fsolve(@new, x0);
% disp(['Omega = ' num2str(x(1)) ' rad'])
% disp(['  Tau = ' num2str(x(2)) ' hrs'])


%% Solve

% 2D system:
delays = [2];
history = [1e-6 0];
tspan = [0, 95];
sol = dde23(@derivs2D, delays, history, tspan);

plot2D(sol,history, 0,0)

function roots = new(x)
    %initialize model parameters:
    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = 0.5;
    Nm = 1.2; 
    k = 0.9;
    Ae = k*de / (d*ke);

    omega = x(1);
    tau = x(2);
    roots = zeros(2,1);

    roots(1) = k*da*cos(omega*tau) + va*da*cos(omega*tau)^2 - va*da*sin(omega*tau)^2-va*omega*sin(omega*tau) - omega*da*sin(omega*tau) - omega^2;
    roots(2) = -k*da*sin(omega*tau) - omega*k - 2*va*da*cos(omega*tau)*sin(omega*tau) - omega*va*cos(omega*tau) - omega*da*cos(omega*tau);
end


function dydt = derivs2D(t, y, Z)
    N_tau = Z(1,1);

    d = 4e-3;
    ke = 5;

    de = 2;
    va = 4.8e-7;
    da = 0.5;
    Nm = 1.2; % might be e-9
    k = 0.9;

    % Ae = k*de / (d*ke);

    % y = [N E A]
    N = y(1);
    A = y(2);
    E = ke/de*A;

    dNdt = k * N * (1 - N_tau / Nm) - d*E*N;
    dAdt = va*N- da*A;

    dydt = [dNdt;
            dAdt];
end