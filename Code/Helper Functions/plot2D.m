
function plot2D(sol, init_conditions, tau)
    figure; 
    
    time = sol.x;
    y = sol.y;
    N = y(1,:);
    A = y(2,:);

    subplot(2,1,1)
    plot(time, N, 'DisplayName','N(t)');
    hold on
    scatter(0, N(1), 'DisplayName', 'N(0)')
    hold off
    legend;
    title('N(t) vs Time')
    xlabel('Time (hrs)')
    ylabel('N(t) (CFU mL^(-1))')
    subplot(2,1,2)
    plot(time, A, 'DisplayName','A(t)');
    hold on
    scatter(0, A(1), 'DisplayName', 'A(0)')
    hold off
    title('A(t) vs Time')
    xlabel('Time (hrs)')
    ylabel('A(t) (nM)')
    legend;

    sgtitle(['$\tau$ = ' num2str(tau) ' hrs', ...
             ' N(0) = ' num2str(init_conditions(1)), ...
             ' A(0) = ' num2str(init_conditions(2))], 'Interpreter', 'latex')
end 