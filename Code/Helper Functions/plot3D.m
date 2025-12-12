
function plot3D(sol, init_conditions, tau)
    figure;
    
    time = sol.x;
    y = sol.y;
    N = y(1,:);
    E = y(2,:);
    A = y(3,:);

    subplot(3,1,1)
    plot(time, N, 'DisplayName','N(t)');
    legend;
    title('N(t) vs Time')
    xlabel('Time (hrs)')
    ylabel('N(t) (CFU mL^(-1))')
    
    subplot(3,1,2)
    plot(time, E, 'DisplayName','E(t)');
    title('E(t) vs Time')
    xlabel('Time (hrs)')
    ylabel('E(t) (nM)')
    legend;

    subplot(3,1,3)
    plot(time, A, 'DisplayName','A(t)');
    title('A(t) vs Time')
    xlabel('Time (hrs)')
    ylabel('A(t) (nM)')
    legend;


    sgtitle(['$\tau$ = ' num2str(tau) ' hrs',... 
             ' N(0) = ' num2str(init_conditions(1)), ' E(0) = ' num2str(init_conditions(2)) ...
             ' A(0) = ' num2str(init_conditions(3))], 'Interpreter', 'latex')
end 