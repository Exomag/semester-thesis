function generatePlotsCaseStudy()
    % Parameters
    params = loadParameters();
    SAMPLING_TIME = params.SAMPLING_TIME;
    HORIZON = params.HORIZON;
    CAR_WIDTH = params.CAR_WIDTH;
    CAR_LENGTH = params.CAR_LENGTH;
    LANE_WIDTH = params.LANE_WIDTH;
    EGO_X0 = params.EGO_X0;
    EGO_Y0 = params.EGO_Y0;
    OBSTACLE_X = params.OBSTACLE_X;
    OBSTACLE_Y = params.OBSTACLE_Y;
    A = params.A;
    DIST_STD = params.DIST_STD;
    DIST_VAR = params.DIST_VAR;

    % Data
    N_trials = 100;
    costs = zeros(N_trials, 2);
    emp_viol = zeros(N_trials, 2);
    t_solver = zeros(N_trials, 2);
    problems = zeros(N_trials, 2);

    % Trials
    for trial = 1:N_trials
        % Moments Robust Approach
        rng(trial)
        params = loadParameters();
        params.N_s = 1000;
        [solData, solution, objective] = solveProblemMRA(params);
        Pviol_joint = calcMonteCarloViolationProbability(params, solData);
        costs(trial, 1) = value(objective);
        emp_viol(trial, 1) = Pviol_joint;
        t_solver(trial, 1) = solution.solvertime;
        problems(trial, 1) = solution.problem;

        % Scenario Approach
        rng(trial);
        params = loadParameters();
        [solData, solution, objective] = solveProblemSA(params);
        Pviol_joint = calcMonteCarloViolationProbability(params, solData);
        costs(trial, 2) = value(objective);
        emp_viol(trial, 2) = Pviol_joint;
        t_solver(trial, 2) = solution.solvertime;
        problems(trial, 2) = solution.problem;
    end

    % Figure - Costs
    costs_boxplot = [costs(problems(:, 1) == 0, 1); costs(problems(:, 2) == 0, 2)];
    groups_boxplot = [zeros(sum(problems(:, 1) == 0), 1); ones(sum(problems(:, 2) == 0), 1)];
    figure();
    hold on;
    grid
    boxplot(costs_boxplot, groups_boxplot, 'Labels', {'MRA', 'SA'}, 'Whisker', Inf)
    ylabel('cost', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    legend('hide')
    save2tikz('plots/PassObs_CostBoxplot')

    % Figure - Empirical violations
    viols_boxplot = [emp_viol(problems(:, 1) == 0, 1); emp_viol(problems(:, 2) == 0, 2)];
    figure();
    hold on;
    grid
    boxplot(viols_boxplot, groups_boxplot, 'Labels', {'MRA', 'SA'}, 'Whisker', Inf)
    ylabel('violation probability', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    legend('hide')
    save2tikz('plots/PassObs_ViolsBoxplot')

    % Figure - Solver time
    time_boxplot = [t_solver(problems(:, 1) == 0, 1); t_solver(problems(:, 2) == 0, 2)];
    figure();
    hold on;
    grid
    boxplot(time_boxplot, groups_boxplot, 'Labels', {'MRA', 'SA'}, 'Whisker', Inf)
    ylabel('solver time (\si{\second})', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    legend('hide')
    save2tikz('plots/PassObs_RuntimeBoxplot')

    % Figure - Initial state
    params = loadParameters();
    paramOpt = calcOptimizationParameters(params);
    x_adv_nom = paramOpt.x_adv_nom;
    figure();
    hold on;
    grid
    plot(x_adv_nom(1, :), x_adv_nom(3, :), 'r-x', 'LineWidth', 1)
    rectangle('Position', [EGO_X0 - CAR_LENGTH / 2, EGO_Y0 - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'b')
    rectangle('Position', [x_adv_nom(1, 1) - CAR_LENGTH / 2, x_adv_nom(3, 1) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'r')
    rectangle('Position', [OBSTACLE_X - CAR_LENGTH / 2, OBSTACLE_Y - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'k')
    xlim_vec = [EGO_X0 - CAR_LENGTH, x_adv_nom(1, 1) + CAR_LENGTH];
    plot(xlim_vec, LANE_WIDTH*ones(2, 1), '--k', 'LineWidth', 1)
    plot(xlim_vec, 0*ones(2, 1), '-k', 'LineWidth', 2)
    plot(xlim_vec, 2*LANE_WIDTH*ones(2, 1), '-k', 'LineWidth', 2)
    xlim(xlim_vec)
    ylim([0, 2 * LANE_WIDTH])
    xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex')
    ylabel('$$x_3$$ (\si{\meter})', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    legend('hide')
    save2tikz('plots/PassObs_InitialState')

    % Figure - Speed uncertainty
    t_axis = SAMPLING_TIME * (0:HORIZON);
    figure();
    hold on;
    grid
    h1 = plot(t_axis, x_adv_nom(2, :), 'b-', 'LineWidth', 1);
    i = 1;
    h2 = plot([t_axis(i), t_axis(i)], [x_adv_nom(2, i) - DIST_STD, x_adv_nom(2, i) + DIST_STD], 'k-+', 'LineWidth', 1);
    for i = 2:length(t_axis)
        plot([t_axis(i), t_axis(i)], [x_adv_nom(2, i) - DIST_STD, x_adv_nom(2, i) + DIST_STD], 'k-+', 'LineWidth', 1)
    end
    legend([h1, h2], {'mean', '$$\pm 1$$ std'}, 'interpreter', 'latex')
    ylim([-26, -14])
    xlabel('time $$t$$ (\si{\second})', 'interpreter', 'latex')
    ylabel('velocity $$\chi_2$$ (\si{\meter\per\second})', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/PassObs_VelDist')

    % Figure - Position uncertainty
    P = diag([0, DIST_VAR]);
    A_x = A(1:2, 1:2);
    x_std = zeros(size(t_axis));
    x_std(1) = 0;
    for i = 2:length(t_axis)
        P = A_x * P * A_x';
        x_std(i+1) = sqrt(P(1, 1));
    end
    figure();
    hold on;
    grid
    h1 = plot(t_axis, x_adv_nom(1, :), 'b-', 'LineWidth', 1);
    i = 1;
    h2 = plot([t_axis(i), t_axis(i)], [x_adv_nom(1, i) - x_std(i), x_adv_nom(1, i) + x_std(i)], 'k-+', 'LineWidth', 1);
    for i = 2:length(t_axis)
        plot([t_axis(i), t_axis(i)], [x_adv_nom(1, i) - x_std(i), x_adv_nom(1, i) + x_std(i)], 'k-+', 'LineWidth', 1)
    end
    legend([h1, h2], {'mean', '$$\pm 1$$ std'}, 'interpreter', 'latex')
    hold off
    xlabel('time $$t$$ (\si{\second})', 'interpreter', 'latex')
    ylabel('position $$\chi_1$$ (\si{\meter})', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/PassObs_PosDist')

    % Figure - MRA (std = 5)
    rng(0)
    params = loadParameters();
    params.N_s = 1000;
    solData = solveProblemMRA(params);
    if isempty(solData)
        error('Infeasible MRA problem.');
    end
    x_val = solData.x_val;
    ind_frames = 1:4:HORIZON + 1;
    N_frames = length(ind_frames);
    figure()
    for i = 1:N_frames
        ind = ind_frames(i);
        subplot(N_frames, 1, i);
        hold on;
        grid
        plot(x_val(1, 1:ind), x_val(3, 1:ind), 'b-x', 'LineWidth', 1)
        plot(x_adv_nom(1, 1:ind), x_adv_nom(3, 1:ind), 'r-x', 'LineWidth', 1)
        rectangle('Position', [x_val(1, ind) - CAR_LENGTH / 2, x_val(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'b')
        rectangle('Position', [x_adv_nom(1, ind) - CAR_LENGTH / 2, x_adv_nom(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'r')
        rectangle('Position', [OBSTACLE_X - CAR_LENGTH / 2, OBSTACLE_Y - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'k')
        xlim_vec = [EGO_X0 - CAR_LENGTH, max(x_val(1, end), x_adv_nom(1, 1)) + CAR_LENGTH];
        plot(xlim_vec, LANE_WIDTH*ones(2, 1), '--k', 'LineWidth', 1)
        plot(xlim_vec, 0*ones(2, 1), '-k', 'LineWidth', 2)
        plot(xlim_vec, 2*LANE_WIDTH*ones(2, 1), '-k', 'LineWidth', 2)
        xlim(xlim_vec)
        ylim([0, 2 * LANE_WIDTH])
        ylabel('$$x_3$$ (\si{\meter})', 'interpreter', 'latex')
        legend('hide')
        set(gca, 'TickLabelInterpreter', 'latex')
    end
    xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex')
    save2tikz('plots/PassObs_FramesMRA')

    % Figure - SA (std = 5)
    rng(0)
    params = loadParameters();
    solData = solveProblemSA(params);
    if isempty(solData)
        error('Infeasible SA problem.');
    end
    x_val = solData.x_val;
    ind_frames = 1:4:HORIZON + 1;
    N_frames = length(ind_frames);
    figure()
    for i = 1:N_frames
        ind = ind_frames(i);
        subplot(N_frames, 1, i);
        hold on;
        grid
        plot(x_val(1, 1:ind), x_val(3, 1:ind), 'b-x', 'LineWidth', 1)
        plot(x_adv_nom(1, 1:ind), x_adv_nom(3, 1:ind), 'r-x', 'LineWidth', 1)
        rectangle('Position', [x_val(1, ind) - CAR_LENGTH / 2, x_val(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'b')
        rectangle('Position', [x_adv_nom(1, ind) - CAR_LENGTH / 2, x_adv_nom(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'r')
        rectangle('Position', [OBSTACLE_X - CAR_LENGTH / 2, OBSTACLE_Y - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'k')
        xlim_vec = [EGO_X0 - CAR_LENGTH, max(x_val(1, end), x_adv_nom(1, 1)) + CAR_LENGTH];
        plot(xlim_vec, LANE_WIDTH*ones(2, 1), '--k', 'LineWidth', 1)
        plot(xlim_vec, 0*ones(2, 1), '-k', 'LineWidth', 2)
        plot(xlim_vec, 2*LANE_WIDTH*ones(2, 1), '-k', 'LineWidth', 2)
        xlim(xlim_vec)
        ylim([0, 2 * LANE_WIDTH])
        ylabel('$$x_3$$ (\si{\meter})', 'interpreter', 'latex')
        legend('hide')
        set(gca, 'TickLabelInterpreter', 'latex')
    end
    xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex')
    save2tikz('plots/PassObs_FramesSA')

    % Figure - MRA (std = 10)
    rng(0)
    params = loadParameters();
    params.DIST_STD = 10;
    params.N_s = 1000;
    solData = solveProblemMRA(params);
    if isempty(solData)
        error('Infeasible MRA problem.');
    end
    x_val = solData.x_val;
    ind_frames = 1:4:HORIZON + 1;
    N_frames = length(ind_frames);
    figure()
    for i = 1:N_frames
        ind = ind_frames(i);
        subplot(N_frames, 1, i);
        hold on;
        grid
        plot(x_val(1, 1:ind), x_val(3, 1:ind), 'b-x', 'LineWidth', 1)
        plot(x_adv_nom(1, 1:ind), x_adv_nom(3, 1:ind), 'r-x', 'LineWidth', 1)
        rectangle('Position', [x_val(1, ind) - CAR_LENGTH / 2, x_val(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'b')
        rectangle('Position', [x_adv_nom(1, ind) - CAR_LENGTH / 2, x_adv_nom(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'r')
        rectangle('Position', [OBSTACLE_X - CAR_LENGTH / 2, OBSTACLE_Y - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'k')
        xlim_vec = [EGO_X0 - CAR_LENGTH, max(x_val(1, end), x_adv_nom(1, 1)) + CAR_LENGTH];
        plot(xlim_vec, LANE_WIDTH*ones(2, 1), '--k', 'LineWidth', 1)
        plot(xlim_vec, 0*ones(2, 1), '-k', 'LineWidth', 2)
        plot(xlim_vec, 2*LANE_WIDTH*ones(2, 1), '-k', 'LineWidth', 2)
        xlim(xlim_vec)
        ylim([0, 2 * LANE_WIDTH])
        ylabel('$$x_3$$ (\si{\meter})', 'interpreter', 'latex')
        legend('hide')
        set(gca, 'TickLabelInterpreter', 'latex')
    end
    xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex')
    save2tikz('plots/PassObs_FramesMRA_2')

    % Figure - SA (std = 10)
    rng(0)
    params = loadParameters();
    params.DIST_STD = 10;
    solData = solveProblemSA(params);
    if isempty(solData)
        error('Infeasible SA problem.');
    end
    x_val = solData.x_val;
    ind_frames = 1:4:HORIZON + 1;
    N_frames = length(ind_frames);
    figure()
    for i = 1:N_frames
        ind = ind_frames(i);
        subplot(N_frames, 1, i);
        hold on;
        grid
        plot(x_val(1, 1:ind), x_val(3, 1:ind), 'b-x', 'LineWidth', 1)
        plot(x_adv_nom(1, 1:ind), x_adv_nom(3, 1:ind), 'r-x', 'LineWidth', 1)
        rectangle('Position', [x_val(1, ind) - CAR_LENGTH / 2, x_val(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'b')
        rectangle('Position', [x_adv_nom(1, ind) - CAR_LENGTH / 2, x_adv_nom(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'r')
        rectangle('Position', [OBSTACLE_X - CAR_LENGTH / 2, OBSTACLE_Y - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'k')
        xlim_vec = [EGO_X0 - CAR_LENGTH, max(x_val(1, end), x_adv_nom(1, 1)) + CAR_LENGTH];
        plot(xlim_vec, LANE_WIDTH*ones(2, 1), '--k', 'LineWidth', 1)
        plot(xlim_vec, 0*ones(2, 1), '-k', 'LineWidth', 2)
        plot(xlim_vec, 2*LANE_WIDTH*ones(2, 1), '-k', 'LineWidth', 2)
        xlim(xlim_vec)
        ylim([0, 2 * LANE_WIDTH])
        ylabel('$$x_3$$ (\si{\meter})', 'interpreter', 'latex')
        legend('hide')
        set(gca, 'TickLabelInterpreter', 'latex')
    end
    xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex')
    save2tikz('plots/PassObs_FramesSA_2')
end
