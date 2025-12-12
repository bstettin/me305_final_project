clear
clc
close all

d = 4e-3;
ke = 5;
de = 2;
va = 4.8e-7;
da = mean([0.274, 1.17]);%0.5;
Nm = mean([1.17, 1.25]) *1e9; %1.2; 
k = mean([0.885 0.936]); %0.9;

%calculate equilibrium points
Ne = k / (k/Nm + d*ke/de*va/da);
Ae = va/da*Ne;
Ee = ke*Ae/de;

% 2D system:
delays = [2.4095];
history = [Ne Ee, 80];
tspan = [0, 18];
sol = dde23(@derivs3D, delays, history, tspan);

plot3D(sol,history, delays(1))

figure;
plot(sol.y(2,:), sol.y(1,:))


function dydt = derivs3D(t, y, Z)
    N_tau = Z(1,1);
    
    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = mean([0.274, 1.17]);%0.5;
    Nm = mean([1.17, 1.25]) *1e9; %1.2; 
    k = mean([0.885 0.936]); %0.9;

    % Ne = k / (k/Nm + d*ke/de*va/da);
    % Ae = va/da*Ne;

    % Ae = k*de / (d*ke);

    % y = [N E A]
    N = y(1);
    E = y(2);
    A = y(3);

    dNdt = k * N * ( 1 -  N / Nm) - d*E*N_tau;
    dEdt = ke*A - de*E;
    dAdt = va*N - da*A;

    dydt = [dNdt;
            dEdt;
            dAdt];
end