clear
clc
close all


%% Determine the roots of the system for the 2D case...
d = 4e-3;
ke = 5;
de = 2;
va = 4.8e-7;
da = 0.5;
Nm = 1.2; 
k = 0.9;

Ae = k*de / (d*ke);
Ee = ke*Ae / de;


x0 = [pi/10, 250];
x = fsolve(@roots, x0);
disp(['Omega = ' num2str(x(1)) ' rad'])
disp(['  Tau = ' num2str(x(2)) ' hrs'])
%% Create resulting plots 

% 2D system:
delays = x(2);
history = [1e-6 1e-3];
tspan = [0, 200];
sol = dde23(@derivs2D, delays, history, tspan);

plot2D(sol,history, x(1),x(2))

%% 2D Roots:
Ae = k*de / (d*ke);
Ne = da/va*Ae;
Ne

%%
%3D system:
delays = x(2);
history = [1e-6 Ee, Ae];
tspan = [0, 200];
sol = dde23(@derivs3D, delays, history, tspan);
plot3D(sol, history, x(1), x(2))



%% Functions to calculate derivatives and evaluate the roots

function dydt = derivs2D(t, y, Z)
    delay1 = Z(:,1);

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

    dNdt = k*N* - d*E*delay1(1);
    dAdt = va*N - da*A;

    dydt = [dNdt;
            dAdt];
end

function dydt = derivs3D(t, y, Z)
    delay1 = Z(:,1);

    d = 4e-3;
    ke = 5;

    de = 2;

    %need to change this to induce oscillations in the 3D system by 
    va = 4.8e-7; 
    da = 0.35;
    
    Nm = 1.2; % might be e-9
    k = 0.9;

    % Ae = k*de / (d*ke);

    % y = [N E A]
    N = y(1);
    E = y(2);
    A = y(3);

    dNdt = k*N* - d*E*delay1(1);
    dAdt = va*N - da*A;
    dEdt = ke*A - de*E;

    dydt = [dNdt;
            dEdt;
            dAdt];

    % dydt = clip(dydt, 0, 1e5);
end





function roots = system2D(x)
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

    roots(1) = k*da - d*ke/de*Ae*cos(omega*tau) + d*ke/de*Ae*omega*sin(omega*tau) + omega^2;
    roots(2) = k*omega - d*ke/de*Ae*sin(omega*tau)-d*ke/de*Ae*omega*cos(omega*tau) - da*omega;
end
