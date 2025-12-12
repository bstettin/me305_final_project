clear
clc
close all

    tspan = [0 400];

Nm = 1.2; % might be e-9

trials = 3;
time = cell(1,trials);
N = cell(1,trials);
E = cell(1,trials);
A = cell(1,trials);
multiplier = [1 1/60 1/3600];%[1 1/60 1/3600];
time_delay = 2; %3.3 seems to be max value can have in here before system breaks

d = 4e-3;
ke = 5;

de = 2;
va = 4.8e-7; %4.8e-7; %4.8e-7; %0.5 induces delay
da = 0.5;
Nm = 1.2; % might be e-9
k = 0.9;

Ee = k/d;
Ne = Nm*da*de*k / (Nm*va*ke*d + da*de*k);
Ae = va*Ne / da;

% Ae, Ne, Ee

for i = 1:length(multiplier)
    delays = time_delay * multiplier(i); %1 second delay ;
    history = [1 0 0];
    sol = dde23(@derivs, delays, history, tspan);
    
    time{i} = num2cell(sol.x);
    y = sol.y;
    N{i} = num2cell(y(1,:));
    E{i} = num2cell(y(2,:));
    A{i} = num2cell(y(3,:));
end 

 
figure;

subplot(3,1,1)
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

subplot(3,1,2)

for i=1:trials
    hold on
    delays = time_delay * multiplier(i);
    t = cell2mat(time{i});
    e = cell2mat(E{i});
    plot(t, e, 'DisplayName', ['$\tau$ = ' num2str(round(delays(1), 4)) ' hours'])
    ylabel('E')
    title('E vs t')
    xlabel('Time (hours)')
    legend('Interpreter','latex')
end

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

    
    % y = [N E A]
    N = y(1);
    E = y(2);
    A = y(3);

    dNdt = k*N*(1 - delay1(1) / Nm) - d*E*N;
    dEdt = ke*A - de*E;
    dAdt = va*N - da*A;

    dydt = [dNdt;
            dEdt;
            dAdt];
end
