clear
clc
close all

    tspan = [0 40];

Nm = 1.2; % might be e-9

trials = 3;
time = cell(1,trials);
N = cell(1,trials);
E = cell(1,trials);
A = cell(1,trials);
multiplier = [1 1/60 1/3600];%[1 1/60 1/3600];
time_delay = [1]; %3.3 seems to be max value can have in here before system breaks

d = 4e-3;
ke = 5;

de = 2;
va = 5;%4.8e-7; %4.8e-7; %4.8e-7; %0.5 induces delay
da = 0.5;
Nm = 1.2; % might be e-9
k = 0.9;

Ee = k/d;
Ne = Nm*da*de*k / (Nm*va*ke*d + da*de*k);
Ae = va*Ne / da;

% Ae, Ne, Ee

for i = 1:length(multiplier)
    delays = time_delay * multiplier(i); %1 second delay ;
    history = [1e-5 100];
    sol = dde23(@derivs, delays, history, tspan);
    
    time{i} = num2cell(sol.x);
    y = sol.y;
    N{i} = num2cell(y(1,:));
    A{i} = num2cell(y(2,:));
end 

 
figure;

subplot(2,1,1)
for i=1:trials
    
    hold on
    delays = time_delay * multiplier(i);
    t = cell2mat(time{i});
    n = cell2mat(N{i});
    
    plot(t, n, 'DisplayName', ['$\tau$ = ' num2str(round(delays(1),4)) ' hours'])
    
end 

yline(Nm, 'DisplayName','Carrying Capacity', 'Color','k', 'LineStyle','--')
ylabel('N')
title('N vs t')
xlabel('Time (hours)')
legend('Interpreter','latex')
hold off

subplot(2,1,2)

subplot(3,1,3)
for i=1:trials
    hold on
    delays = time_delay * multiplier(i);
    t = cell2mat(time{i});
    a = cell2mat(A{i});
    plot(t, a, 'DisplayName', ['$\tau$ = ' num2str(round(delays(1), 4)) ' hours'])
    ylabel('A')
    title('A vs t')
    xlabel('Time (hours)')
    legend('Interpreter','latex')
end
% 
% figure;
% plot()



function dydt = derivs(t, y, Z)
    delay1 = Z(:,1);

    d = 4e-3;
    ke = 5;

    

    de = 2;
    va = 4.8e-7;%4.8e-7; %4.8e-7; %4.8e-7; %0.5 induces delay
    da = 0.5;
    Nm = 1.2; % might be e-9
    k = 0.9;

    
    k = k/3;

    
    % y = [N E A]
    N = y(1);
    A = y(2);

    dNdt = k*delay1(1) - d*ke/de*A*N;
    % dEdt = ke*A - de*E;
    dAdt = va*N - da*A;

    dydt = [dNdt;
            dAdt];
end
