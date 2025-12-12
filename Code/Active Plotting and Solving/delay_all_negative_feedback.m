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
Ne, Ae

x0 = [0.2, 1];
% x0 = [0.25,2]
options = optimoptions('fsolve', ...
    'MaxFunctionEvaluations', 1e5, ... % Increase the function evaluation limit (e.g., to 2000)
    'MaxIterations', 200, ...            % Increase the maximum iterations (e.g., to 500)
    'Display', 'iter');

x = fsolve(@all_negative_roots, x0, options);
disp(['Omega = ' num2str(x(1)) ' rad'])
disp(['  Tau = ' num2str(x(2)) ' hrs'])



%%
%Solve

% 2D system:
% delays = [x(2)];
% history = [Ne 80];
% tspan = [0, 20];
% sol = dde23(@derivs2D, delays, history, tspan);
% 
% plot2D(sol,history, delays(1));
% 
% figure;
% plot(sol.y(2,:), sol.y(1,:))
% xlabel('A(t)');
% ylabel('N(t)');
% title('N(t) vs A(t)')

A_start = 80;
figure;
for i = 0.01:0.5:2
    delays = [i];
    history = @(t)[Ne A_start];

    tend = 15;
    if i < 1.5
        tend = 30;
    end 
    tspan = [0, tend];

    hist = [Ne A_start];
    sol = dde23(@derivs2D, delays, history, tspan);

    time = sol.x;
    y = sol.y;
    N = y(1,:);
    A = y(2,:);

    subplot(2,1,1)
    plot(time, N, 'DisplayName',['\tau =' num2str(i) ' hrs']);
    hold on
    legend('Location','southeast');
    title('N(t) vs Time')
    xlabel('Time (hrs)')
    ylabel('N(t) (CFU mL^{-1})')

    subplot(2,1,2)
    plot(time, A, 'DisplayName',['\tau = ' num2str(i) ' hrs']);
    hold on
    title('A(t) vs Time')
    xlabel('Time (hrs)')
    ylabel('A(t) (nM mL^{-1})')
    legend('Location','southeast');
end
sgtitle('Two State Delay System 2: System Dynamics with Various Time Delays')

function roots = all_negative_roots(x)
    %initialize model parameters:
    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = mean([0.274, 1.17]);%0.5;
    Nm = mean([1.17, 1.25]) *1e9; %1.2; 
    k = mean([0.885 0.936]); %0.9;
    
    Ne = k / (k/Nm + d*ke/de*va/da);
    Ae = va/da*Ne;


    omega = x(1);
    tau = x(2);
    roots = zeros(2,1);

    roots(1) = -k*da + k*Ne/Nm*da*k*da*cos(omega*tau) + k*omega*sin(omega*tau) - omega^2 + va*Ne*d*ke/de;
    roots(2) = -k*da*sin(omega*tau) + omega*da - k*omega + k*omega*Ne/Nm + k*omega*cos(omega*tau);
end


function dydt = derivs2D(t, y, Z)
    N_tau = Z(1,1);

    d = 4e-3;
    ke = 5;
    de = 2;
    va = 4.8e-7;
    da = mean([0.274, 1.17]);%0.5;
    Nm = mean([1.17, 1.25]) *1e9; %1.2; 
    k = mean([0.885 0.936]); %0.9;

    % Ae = k*de / (d*ke);
    % y = [N E A]

    N = y(1);
    A = y(2);
    E = ke/de*A;

    dNdt = k * N - (k*N / Nm + d*E)*N_tau;
    dAdt = va*N- da*A;

    dydt = [dNdt;
            dAdt];
end