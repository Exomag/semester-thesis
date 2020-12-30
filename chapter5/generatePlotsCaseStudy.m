function generatePlotsCaseStudy()
    % Parameters
    params = loadParameters();
    HORIZON = params.HORIZON;
    CAR_WIDTH = params.CAR_WIDTH;
    CAR_LENGTH = params.CAR_LENGTH;
    LANE_WIDTH = params.LANE_WIDTH;

    % Data
    N_trials = 100;
    trialCosts = zeros(N_trials, 1);
    trialViols = zeros(N_trials, 1);
    trialTimes = zeros(N_trials, 1);
    trialFeas = zeros(N_trials, 1);

    % Trials
    demo_mode = false;
    for trial = 1:N_trials
        rng(trial)
        simData = solveProblemRecedingHorizon(params, demo_mode);
        trialCosts(trial, 1) = simData.totalCost;
        trialViols(trial, 1) = double(simData.collision);
        trialTimes(trial, 1) = simData.totalTime;
        trialFeas(trial, 1) = double(simData.feasibleProblem);
    end

    % Demo
    demo_mode = true;
    rng(0)
    simData = solveProblemRecedingHorizon(params, demo_mode);
    if ~simData.feasibleProblem
        error('Infeasible receding horizon problem.');
    end
    tHorizon = simData.tHorizon;
    ego_StatePlan = simData.ego_StatePlan;
    adv1_StateNoNoise = simData.adv1_StateNoNoise;
    ego_State = simData.x_val;
    ego_Input = simData.u_val;
    adv1_State = simData.adv1_NomState;
    adv2_State = simData.adv2_NomState;
    adv2_mean_est = simData.adv2_mean_est;
    adv2_var_est = simData.adv2_var_est;

    % Figure - Costs
    trialCosts_boxplot = trialCosts(trialFeas == 1);
    figure();
    hold on;
    grid
    boxplot(trialCosts_boxplot, 'Whisker', Inf)
    ylabel('cost', 'interpreter', 'latex')
    set(gca, 'XTickLabel', {' '})
    set(gca, 'TickLabelInterpreter', 'latex')
    legend('hide')
    save2tikz('plots/Overtaking_BoxCosts')

    % Figure - Solver time
    trialtimes_boxplot = trialTimes(trialFeas == 1);
    figure();
    hold on;
    grid
    boxplot(trialtimes_boxplot, 'Whisker', Inf)
    ylabel('solver time (\si{\second})', 'interpreter', 'latex')
    set(gca, 'XTickLabel', {' '})
    set(gca, 'TickLabelInterpreter', 'latex')
    legend('hide')
    save2tikz('plots/Overtaking_BoxTimes')

    % Figure - Initial state
    figure();
    hold on;
    grid
    plot(adv1_StateNoNoise(1, :), adv1_StateNoNoise(3, :), 'm-x', 'LineWidth', 1)
    plot(adv2_State(1, :), adv2_State(3, :), 'r-x', 'LineWidth', 1)
    rectangle('Position', [ego_State(1, 1) - CAR_LENGTH / 2, ego_State(3, 1) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'b')
    rectangle('Position', [adv1_StateNoNoise(1, 1) - CAR_LENGTH / 2, adv1_StateNoNoise(3, 1) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'm')
    rectangle('Position', [adv2_State(1, 1) - CAR_LENGTH / 2, adv2_State(3, 1) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'r')
    xlim_vec = [min([min(ego_State(1, :)), min(adv1_StateNoNoise(1, :)), min(adv2_State(1, :))]), max([max(ego_State(1, :)), max(adv1_StateNoNoise(1, :)), max(adv2_State(1, :))])] + [-CAR_LENGTH, CAR_LENGTH];
    plot(xlim_vec, LANE_WIDTH*ones(2, 1), '--k', 'LineWidth', 1)
    plot(xlim_vec, 0*ones(2, 1), '-k', 'LineWidth', 2)
    plot(xlim_vec, 2*LANE_WIDTH*ones(2, 1), '-k', 'LineWidth', 2)
    xlim(xlim_vec)
    ylim([0, 2 * LANE_WIDTH])
    xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex')
    ylabel('$$x_3$$ (\si{\meter})', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    legend('hide')
    save2tikz('plots/Overtaking_InitState')

    % Figure - Position estimate
    figure();
    hold on;
    grid
    h1 = plot(tHorizon, adv2_mean_est(1, :), '-b', 'LineWidth', 1);
    h2 = plot(tHorizon, adv2_State(1, :), '-r', 'LineWidth', 1);
    h3 = plot(tHorizon, adv2_mean_est(1, :)+sqrt(adv2_var_est(1, :)), '--k', 'LineWidth', 1);
    plot(tHorizon, adv2_mean_est(1, :)-sqrt(adv2_var_est(1, :)), '--k', 'LineWidth', 1)
    legend([h1, h2, h3], {'estimate', 'actual', '$$\pm 1$$ std'}, 'interpreter', 'latex')
    xlabel('time $$t$$ (\si{\second})', 'interpreter', 'latex')
    ylabel('position $$\chi_1$$ (\si{\meter})', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/Overtaking_PosEst')

    % Figure - Velocity estimate
    figure();
    hold on;
    grid
    h1 = plot(tHorizon, adv2_mean_est(2, :), '-b', 'LineWidth', 1);
    h2 = plot(tHorizon, adv2_State(2, :), '-r', 'LineWidth', 1);
    h3 = plot(tHorizon, adv2_mean_est(2, :)+sqrt(adv2_var_est(2, :)), '--k', 'LineWidth', 1);
    plot(tHorizon, adv2_mean_est(2, :)-sqrt(adv2_var_est(2, :)), '--k', 'LineWidth', 1)
    legend([h1, h2, h3], {'estimate', 'actual', '$$\pm 1$$ std'}, 'interpreter', 'latex')
    xlabel('time $$t$$ (\si{\second})', 'interpreter', 'latex')
    ylabel('velocity $$\chi_2$$ (\si{\meter\per\second})', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/Overtaking_VelEst')

    % Figure - Speed
    figure();
    hold on;
    grid
    plot(tHorizon, ego_State(2, :), '-b', 'LineWidth', 1)
    plot(tHorizon, ego_State(4, :), '-r', 'LineWidth', 1)
    xlabel('time $$t$$ (\si{\second})', 'interpreter', 'latex')
    ylabel('velocity (\si{\meter\per\second})', 'interpreter', 'latex')
    legend({'$$x_2$$', '$$x_4$$'}, 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/Overtaking_Velocity')

    % Figure - Acceleration
    figure();
    hold on;
    grid
    plot(tHorizon(1:end - 1), ego_Input(1, :), '-b', 'LineWidth', 1)
    plot(tHorizon(1:end - 1), ego_Input(2, :), '-r', 'LineWidth', 1)
    xlabel('time $$t$$ (\si{\second})', 'interpreter', 'latex')
    ylabel('acceleration (\si{\meter\per\square\second})', 'interpreter', 'latex')
    legend({'$$u_1$$', '$$u_2$$'}, 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/Overtaking_Acceleration')

    % Figure - Frames
    ind_frames = 1:4:HORIZON + 1;
    N_frames = length(ind_frames);
    figure()
    for i = 1:N_frames
        ind = ind_frames(i);
        subplot(N_frames, 1, i);
        hold on;
        grid
        plot(ego_State(1, 1:ind), ego_State(3, 1:ind), 'b-x', 'LineWidth', 1)
        plot(ego_StatePlan{ind}(1, :), ego_StatePlan{ind}(3, :), 'g-x', 'LineWidth', 1)
        plot(adv1_State(1, 1:ind), adv1_State(3, 1:ind), 'm-x', 'LineWidth', 1)
        plot(adv2_State(1, 1:ind), adv2_State(3, 1:ind), 'r-x', 'LineWidth', 1)
        rectangle('Position', [ego_State(1, ind) - CAR_LENGTH / 2, ego_State(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'b')
        rectangle('Position', [adv1_State(1, ind) - CAR_LENGTH / 2, adv1_State(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'm')
        rectangle('Position', [adv2_State(1, ind) - CAR_LENGTH / 2, adv2_State(3, ind) - CAR_WIDTH / 2, CAR_LENGTH, CAR_WIDTH], 'FaceColor', 'r')
        xlim_vec = [min([min(ego_State(1, :)), min(adv1_State(1, :)), min(adv2_State(1, :))]), max([max(ego_State(1, :)), max(adv1_State(1, :)), max(adv2_State(1, :))])] + [-CAR_LENGTH, CAR_LENGTH];
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
    save2tikz('plots/Overtaking_Frames')
end
