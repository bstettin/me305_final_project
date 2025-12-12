clear
clc
close all

d = 4e-3;
ke = 5;
de = 2;
va = 4.8e-7;
da = mean([0.274, 1.17]);%0.5;
Nm = mean([1.17, 1.25]) * 1e9; %1.2; 
k = mean([0.885 0.936]); %0.9;

%calculate equilibrium points
Ne = k / (k/Nm + d*ke/de*va/da);
Ae = va/da*Ne;


x0 = [0 1];
x = fsolve(@single_delay, x0);
disp(['Omega = ' num2str(x(1)) ' rad'])
disp(['  Tau = ' num2str(x(2)) ' hrs'])



A_start = 80;
% 2D system:
delays = [x(2)];
history = @(t)[Ne A_start];
tspan = [0, 15];

hist = [Ne A_start];
sol = dde23(@derivs2D, delays, history, tspan);

plot2D(sol,hist, delays(1))

figure;
hold on
plot(sol.y(2,:), sol.y(1,:))
scatter(A_start, Ne)
hold off

% A_start = 80;
% figure;
% for i = 0.01:0.5:2
%     delays = [i];
%     history = @(t)[Ne A_start];
% 
%     tend = 20;
%     if i < 1.5
%         tend = 50;
%     end 
%     tspan = [0, tend];
% 
%     hist = [Ne A_start];
%     sol = dde23(@derivs2D, delays, history, tspan);
% 
%     time = sol.x;
%     y = sol.y;
%     N = y(1,:);
%     A = y(2,:);
% 
%     subplot(2,1,1)
%     plot(time, N, 'DisplayName',['\tau =' num2str(i) ' hrs']);
%     hold on
%     legend('Location','southeast');
%     title('N(t) vs Time')
%     xlabel('Time (hrs)')
%     ylabel('N(t) (CFU mL^(-1))')
% 
%     subplot(2,1,2)
%     plot(time, A, 'DisplayName',['\tau = ' num2str(i) ' hrs']);
%     hold on
%     title('A(t) vs Time')
%     xlabel('Time (hrs)')
%     ylabel('A(t) (nM)')
%     legend('Location','southeast');
% 
% end
% sgtitle('Two State Delay System 1: System Dynamics with Various Time Delays')


function roots = single_delay(x)
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

    roots(1) = -da*k + 2*k*Ne/Nm*da + d*ke/de*Ae*cos(omega*tau)*da - d*ke/de*Ae*sin(omega*tau)*omega-omega^2 + Ne*d*ke/de*va;
    roots(2) = -d*ke/de*Ae*sin(omega*tau)*da + omega*da - k*omega+2*k*Ne/Nm*omega+d*ke/de*Ae*cos(omega*tau);
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

    % Ne = k / (k/Nm + d*ke/de*va/da);
    % Ae = va/da*Ne;

    % Ae = k*de / (d*ke);

    % y = [N E A]

    N = y(1);
    A = y(2);
    E = ke/de*A;

    dNdt = k * N * ( 1 -  N / Nm) - d*E*N_tau;
    dAdt = va*N - da*A;

    dydt = [dNdt;
            dAdt];
end