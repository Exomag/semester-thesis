function generatePlotsExample()
    % Parameters
    N_ax = 1:1:100;
    e_ax = linspace(1e-3, 5e-1, 100);
    N = 20;
    e = 1e-2;

    % Conservativeness values
    conservativeness_N = (1 - e ./ N_ax).^N_ax - (1 - e);
    conservativeness_N = 100 * conservativeness_N / (1 - e);
    conservativeness_e = (1 - e_ax ./ N).^N - (1 - e_ax);
    conservativeness_e = 100 * conservativeness_e ./ (1 - e_ax);

    % Conservativeness as a function of number of samples
    figure();
    hold on;
    grid on;
    plot(N_ax, conservativeness_N, '-k', 'LineWidth', 1)
    ylim([0, 5e-3])
    title('$$\epsilon = 0.01$$', 'interpreter', 'latex')
    xlabel('$$N$$', 'interpreter', 'latex')
    ylabel('conservativeness (\%)', 'interpreter', 'latex')
    legend('hide')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/BooleInequalityExample_N')

    % Conservativeness as a function of number of safety parameter epsilon
    figure();
    hold on;
    grid on;
    plot(e_ax, conservativeness_e, '-k', 'LineWidth', 1)
    ylim([0, 25])
    title('$$N = 20$$', 'interpreter', 'latex')
    xlabel('$$\epsilon$$', 'interpreter', 'latex')
    ylabel('conservativeness (\%)', 'interpreter', 'latex')
    legend('hide')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/BooleInequalityExample_epsilon')
end
