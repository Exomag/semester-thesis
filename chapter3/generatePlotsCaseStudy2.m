function generatePlotsCaseStudy2()
    % Parameters
    params = loadParameters();
    params.dist_type = 'unif';
    Ts = params.Ts;
    c_x = params.c_x;
    c_y = params.c_y;
    l_x = params.l_x;
    N = params.N;

    % Data
    N_trials = 100;
    costs = zeros(N_trials, 2);
    emp_viol = zeros(N_trials, 2);
    t_yalmip = zeros(N_trials, 2);
    t_solver = zeros(N_trials, 2);
    problems = zeros(N_trials, 2);

    % Trials
    for trial = 1:N_trials
        % Moments Robust Approach
        rng(trial)
        [solData, sol, objective] = solveProblemMRA(params);
        if ~isempty(solData)
            Pviol_joint = calcMonteCarloViolationProbability(params, solData);
        else
            Pviol_joint = NaN;
        end
        costs(trial, 1) = value(objective);
        emp_viol(trial, 1) = Pviol_joint;
        t_yalmip(trial, 1) = sol.yalmiptime;
        t_solver(trial, 1) = sol.solvertime;
        problems(trial, 1) = sol.problem;

        % Scenario Approach
        rng(trial)
        [solData, sol, objective] = solveProblemScenario(params);
        if ~isempty(solData)
            Pviol_joint = calcMonteCarloViolationProbability(params, solData);
        else
            Pviol_joint = NaN;
        end
        costs(trial, 2) = value(objective);
        emp_viol(trial, 2) = Pviol_joint;
        t_yalmip(trial, 2) = sol.yalmiptime;
        t_solver(trial, 2) = sol.solvertime;
        problems(trial, 2) = sol.problem;
    end

    % Figure - Costs
    costs_boxplot = [costs(problems(:, 1) == 0, 1); costs(problems(:, 2) == 0, 2)];
    groups_boxplot = [zeros(sum(problems(:, 1) == 0), 1); ones(sum(problems(:, 2) == 0), 1)];
    figure()
    boxplot(costs_boxplot, groups_boxplot, 'Labels', {'MRA', 'SA'}, 'Whisker', Inf)
    grid
    ylabel('cost', 'interpreter', 'latex')
    legend('hide')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/jointCC_costs_boxplot_case2')

    % Figure - Empirical violations
    viols_boxplot = [emp_viol(problems(:, 1) == 0, 1); emp_viol(problems(:, 2) == 0, 2)];
    figure()
    epsilon = params.epsilon;
    boxplot(viols_boxplot, groups_boxplot, 'Labels', {'MRA', 'SA'}, 'Whisker', Inf)
    hold on
    xl = xlim;
    plot([xl(1), xl(end)], [epsilon, epsilon], 'k--')
    hold off
    ylim([-0.001, epsilon * 1.05])
    grid
    legend('\epsilon', 'interpreter', 'latex')
    ylabel('violation probability', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/jointCC_viols_boxplot_case2')

    % Figure - Solver time
    time_boxplot = [t_solver(problems(:, 1) == 0, 1); t_solver(problems(:, 2) == 0, 2)];
    figure()
    boxplot(time_boxplot, groups_boxplot, 'Labels', {'MRA', 'SA'}, 'Whisker', Inf)
    grid
    ylabel('solver time (\si{\second})', 'interpreter', 'latex')
    legend('hide')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/jointCC_times_boxplot_case2')

    % Figure - Nonmoving frame (SA)
    rng(0)
    [solData, ~, ~, samples, k_active] = solveProblemScenario(params);
    if isempty(solData)
        error('Infeasible SA problem.');
    end
    x1_abs = solData.x1_abs;
    y1_abs = solData.y1_abs;
    x2_abs = solData.x2_abs;
    y2_abs = solData.y2_abs;
    samples_lim = zeros(N, 2);
    for i = 1:N
        samples_lim(i, :) = [min(samples{i}), max(samples{i})];
    end
    samples_lim = [0, 0; samples_lim];
    figure()
    plot(x1_abs, x2_abs, '-bx', 'LineWidth', 1)
    hold on
    plot(y1_abs, y2_abs, '-rx', 'LineWidth', 1)
    plot(y1_abs+samples_lim(:, 2), y2_abs, '--r', 'LineWidth', 1)
    plot([-l_x, -l_x], [min(y2_abs), max(y2_abs)], '-k', 'LineWidth', 2)
    plot([l_x, l_x], [min(y2_abs), max(y2_abs)], '-k', 'LineWidth', 2)
    plot([0, 0], [min(y2_abs), max(y2_abs)], '--k', 'LineWidth', 1)
    hold off
    xlim([-l_x, l_x])
    ylim([min(y2_abs), max(y2_abs)])
    grid
    xlabel('$$x_1$$ (\si{{\metre}})', 'interpreter', 'latex')
    ylabel('$$x_3$$ (\si{{\metre}})', 'interpreter', 'latex')
    legend({'ego car', 'adversary car', 'maximum uncertainty'}, 'Location', 'northwest', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/jointCC_SA_nonmoving_case2')

    % Figure - Moving frame (SA)
    figure()
    y2 = solData.y2;
    indDraw = [0; k_active(1); k_active(round(length(k_active) / 2 + 1)); k_active(end) + 1; N];
    nDraw = length(indDraw);
    for i = 1:nDraw
        iDraw = indDraw(i) + 1;
        subplot(1, nDraw, i)
        hold on
        plot(x1_abs(1:iDraw), zeros(iDraw, 1), '-bx', 'LineWidth', 1)
        rectangle('Position', [x1_abs(iDraw) - c_x / 2, -c_y / 2, c_x, c_y], 'LineWidth', 1)
        plot(y1_abs(1:iDraw), y2(1:iDraw), '-rx', 'LineWidth', 1)
        rectangle('Position', [y1_abs(iDraw) - c_x / 2, y2(iDraw) - c_y / 2, c_x, c_y], 'LineWidth', 1)
        plot([-l_x, -l_x], [min(y2), max(y2)], '-k', 'LineWidth', 2)
        plot([l_x, l_x], [min(y2), max(y2)], '-k', 'LineWidth', 2)
        plot([0, 0], [min(y2), max(y2)], '--k', 'LineWidth', 1)
        hold off
        title(['t = ', num2str(Ts * (iDraw - 1)), 's'])
        xlim([-l_x, l_x])
        ylim([min(y2), max(y2)])
        grid
        title(['$$t = \SI{', num2str(Ts * (iDraw - 1)), '}{\second}$$'], 'interpreter', 'latex')
        xlabel('$$x_1$$ (\si{{\metre}})', 'interpreter', 'latex')
        if i == 1
            ylabel('$$x_3$$ (\si{{\metre}})', 'interpreter', 'latex')
        end
        legend('hide')
        set(gca, 'TickLabelInterpreter', 'latex')
    end
    save2tikz('plots/jointCC_SA_moving_case2')

    % Figure - Nonmoving frame (MRA)
    rng(0)
    [solData, ~, ~, samples, k_active] = solveProblemScenario(params);
    if isempty(solData)
        error('Infeasible SA problem.');
    end
    x1_abs = solData.x1_abs;
    y1_abs = solData.y1_abs;
    x2_abs = solData.x2_abs;
    y2_abs = solData.y2_abs;
    samples_lim = zeros(N, 2);
    for i = 1:N
        samples_lim(i, :) = [min(samples{i}), max(samples{i})];
    end
    samples_lim = [0, 0; samples_lim];
    figure()
    plot(x1_abs, x2_abs, '-bx', 'LineWidth', 1)
    hold on
    plot(y1_abs, y2_abs, '-rx', 'LineWidth', 1)
    plot(y1_abs+[0; ones(length(y1_abs) - 1, 1) * samples_lim(2, 2)], y2_abs, '--r', 'LineWidth', 1)
    plot([-l_x, -l_x], [min(y2_abs), max(y2_abs)], '-k', 'LineWidth', 2)
    plot([l_x, l_x], [min(y2_abs), max(y2_abs)], '-k', 'LineWidth', 2)
    plot([0, 0], [min(y2_abs), max(y2_abs)], '--k', 'LineWidth', 1)
    hold off
    xlim([-l_x, l_x])
    ylim([min(y2_abs), max(y2_abs)])
    grid
    xlabel('$$x_1$$ (\si{{\metre}})', 'interpreter', 'latex')
    ylabel('$$x_3$$ (\si{{\metre}})', 'interpreter', 'latex')
    legend({'ego car', 'adversary car', 'maximum uncertainty'}, 'Location', 'northwest', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/jointCC_MRA_nonmoving_case2')

    % Figure - Moving frame (MRA)
    figure()
    y2 = solData.y2;
    indDraw = [0; k_active(1); k_active(round(length(k_active) / 2 + 1)); k_active(end) + 1; N];
    nDraw = length(indDraw);
    for i = 1:nDraw
        iDraw = indDraw(i) + 1;
        subplot(1, nDraw, i)
        hold on
        plot(x1_abs(1:iDraw), zeros(iDraw, 1), '-bx', 'LineWidth', 1)
        rectangle('Position', [x1_abs(iDraw) - c_x / 2, -c_y / 2, c_x, c_y], 'LineWidth', 1)
        plot(y1_abs(1:iDraw), y2(1:iDraw), '-rx', 'LineWidth', 1)
        rectangle('Position', [y1_abs(iDraw) - c_x / 2, y2(iDraw) - c_y / 2, c_x, c_y], 'LineWidth', 1)
        plot([-l_x, -l_x], [min(y2), max(y2)], '-k', 'LineWidth', 2)
        plot([l_x, l_x], [min(y2), max(y2)], '-k', 'LineWidth', 2)
        plot([0, 0], [min(y2), max(y2)], '--k', 'LineWidth', 1)
        hold off
        title(['t = ', num2str(Ts * (iDraw - 1)), 's'])
        xlim([-l_x, l_x])
        ylim([min(y2), max(y2)])
        grid
        title(['$$t = \SI{', num2str(Ts * (iDraw - 1)), '}{\second}$$'], 'interpreter', 'latex')
        xlabel('$$x_1$$ (\si{{\metre}})', 'interpreter', 'latex')
        if i == 1
            ylabel('$$x_3$$ (\si{{\metre}})', 'interpreter', 'latex')
        end
        legend('hide')
        set(gca, 'TickLabelInterpreter', 'latex')
    end
    save2tikz('plots/jointCC_MRA_moving_case2')
end
