function demoUniform()
    % Parameters
    rng(0)
    beta = 10^-3;
    epsilon = 0.1;
    N_trials = 1e2;
    N_MonteCarlo = 1e7;

    % Distribution parameters
    dist_mean = -0.5;
    dist_var = 1 / 12;
    dist_type = 'unif';
    x_star = sqrt(3) * sqrt(betainv(1 - 2 * epsilon, 1 / 2, 1)) * sqrt(dist_var) + dist_mean;

    % Trials
    x_scen = zeros(N_trials, 1);
    x_dist = zeros(N_trials, 1);
    x_dist_known = zeros(N_trials, 1);
    x_cvar = zeros(N_trials, 1);
    for i = 1:N_trials
        rng(i)

        % Scenario approach
        N_scen = ceil(2/epsilon*log(1 / beta));
        samples_scen = generateSamples(N_scen, dist_mean, dist_var, dist_type);
        x_scen(i) = ScenarioApproach(samples_scen);

        % Distributionally Robust approach
        N_dist = N_scen;
        samples_dist = samples_scen(1:N_dist);
        x_dist(i) = DistRobust(samples_dist, epsilon, beta);

        % Moments Robust approach
        N_dist_known = N_scen;
        samples_dist_known = samples_scen(1:N_dist_known);
        x_dist_known(i) = KnownDist(samples_dist_known, epsilon, beta);

        % CVaR approach
        N_cvar = N_scen;
        samples_cvar = samples_scen(1:N_cvar);
        x_cvar(i) = CvarApproach(samples_cvar, epsilon);
    end

    % Calculate violation probability of each solution via Monte Carlo sampling
    rng(0)
    trials_MonteCarlo = generateSamples(N_MonteCarlo, dist_mean, dist_var, dist_type);
    pViol_scen = zeros(N_trials, 1);
    pViol_dist = zeros(N_trials, 1);
    pViol_dist_known = zeros(N_trials, 1);
    pViol_cvar = zeros(N_trials, 1);
    for i = 1:N_trials
        pViol_scen(i) = sum(x_scen(i) < trials_MonteCarlo) / N_MonteCarlo;
        pViol_dist(i) = sum(x_dist(i) < trials_MonteCarlo) / N_MonteCarlo;
        pViol_dist_known(i) = sum(x_dist_known(i) < trials_MonteCarlo) / N_MonteCarlo;
        pViol_cvar(i) = sum(x_cvar(i) < trials_MonteCarlo) / N_MonteCarlo;
    end

    % Optimization results boxplot 1
    figure();
    hold on;
    grid
    boxplot([x_scen, x_dist_known, x_cvar], 'Whisker', Inf, 'Labels', {'SA', 'MRA', 'CVaRA'})
    xl = xlim;
    plot([xl(1), xl(end)], [x_star, x_star], 'k--')
    ylim([x_star * 1.05, 0])
    ylabel('optimization result', 'interpreter', 'latex')
    legend({'$$x^*$$'}, 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/CCProg_Unif_SolBoxplot1')

    % Optimization results boxplot 2
    figure();
    hold on;
    grid
    boxplot(x_dist, 'Whisker', Inf, 'Labels', {'DRA'})
    ylabel('optimazation result', 'interpreter', 'latex')
    legend('hide')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/CCProg_Unif_SolBoxplot2')

    % Violation probabilities boxplot
    figure();
    hold on;
    grid
    boxplot([pViol_scen, pViol_dist, pViol_dist_known, pViol_cvar], 'Whisker', Inf, 'Labels', {'SA', 'DRA', 'MRA', 'CVaRA'})
    xl = xlim;
    plot([xl(1), xl(end)], [epsilon, epsilon], 'k--')
    hold off
    ylim([0, epsilon * 1.05])
    ylabel('violation probability', 'interpreter', 'latex')
    legend({'$$\epsilon$$'}, 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/CCProg_Unif_ViolBoxplot')

    % Comparison plot
    figure();
    rng(0)
    R = 1;
    r1 = (R / sqrt(N_scen)) * (2 + sqrt(2 * log(2 / beta)));
    r1_alt = abs(dist_mean*(1 / (beta^(1 / N_scen)) - 1));
    hold on;
    grid
    xl = [-1.05, 0.05];
    plot(xl(1):0.001:xl(2), unifpdf(xl(1):0.001:xl(2), -1, 0), 'k-', 'LineWidth', 1)
    y_min = -0.3;
    y_max = 1.2;
    plot([mean(x_star), mean(x_star)], [y_min, y_max], '-*', 'LineWidth', 1)
    plot([mean(x_scen), mean(x_scen)], [y_min, y_max], '-*', 'LineWidth', 1)
    plot([mean(x_dist_known), mean(x_dist_known)], [y_min, y_max], '-*', 'LineWidth', 1)
    plot([mean(x_cvar), mean(x_cvar)], [y_min, y_max], '-*', 'LineWidth', 1)
    plot([dist_mean - r1_alt, dist_mean + r1_alt], [-0.2, -0.2], '-+', 'LineWidth', 1)
    plot([dist_mean - r1, dist_mean + r1], [-0.1, -0.1], '-+', 'LineWidth', 1)
    xlim(xl)
    ylim([y_min, y_max])
    xlabel('$$x$$', 'interpreter', 'latex')
    legend({'$$\mathrm{Unif}[-1,0]$$', '$$x^*$$', 'SA', 'MRA', 'CVaRA', '$$r_1$$', '$$t_1$$'}, 'Location', 'NorthWest', 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/CCProg_Unif_Comparison')
end

function x = ScenarioApproach(samples)
    f = 1;
    A = -ones(size(samples));
    b = -samples;
    options = optimoptions('linprog', 'Display', 'off');
    x = linprog(f, A, b, [], [], [], [], options);
end

function x = DistRobust(samples, epsilon, beta)
    N = length(samples);
    mean_est = mean(samples);
    var_est = var(samples);
    R = 1;
    r1 = (R / sqrt(N)) * (2 + sqrt(2 * log(2 / beta)));
    r2 = (2 * R^2 / sqrt(N)) * (2 + sqrt(2 * log(4 / beta)));
    fun = @(x) x;
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e3, 'Display', 'off');
    x = fmincon(fun, 0, [], [], [], [], [], [], @(x)DistRobustConstr(x, epsilon, mean_est, var_est, r1, r2), options);
end

function x = KnownDist(samples, epsilon, beta)
    beta = beta / 2;
    N = length(samples);
    mean_est = (min(samples) + max(samples)) / 2;
    var_est = (max(samples) - min(samples))^2 / 12;
    r1 = abs(mean_est*(1 / (beta^(1 / N)) - 1));
    r2 = abs(var_est*((1 / (beta^(1 / N)))^2 - 1));
    fun = @(x) x;
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e3, 'Display', 'off');
    x = fmincon(fun, 0, [], [], [], [], [], [], @(x)KnownDistConstr(x, epsilon, mean_est, var_est, r1, r2), options);
end

function x = CvarApproach(samples, epsilon)
    fun = @(x) x;
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e3, 'Display', 'off');
    x = fmincon(fun, 0, [], [], [], [], [], [], @(x)CVarConstr(x, epsilon, samples), options);
end

function [c, ceq] = DistRobustConstr(x, epsilon, mean_est, var_est, r1, r2)
    c = sqrt((1-epsilon)/epsilon) * sqrt(var_est+r2) - x + mean_est + r1;
    ceq = [];
end

function [c, ceq] = KnownDistConstr(x, epsilon, mean_est, var_est, r1, r2)
    c = sqrt(3) * sqrt(betainv(1 - 2 * epsilon, 1 / 2, 1)) * sqrt(var_est+r2) - x + mean_est + r1;
    ceq = [];
end

function [c, ceq] = CVarConstr(x, epsilon, samples)
    alpha0 = 0;
    g_vec = samples - x;
    [~, c] = fminsearch(@(alpha)CVarLinSearch(alpha, epsilon, g_vec), alpha0);
    ceq = [];
end

function y = CVarLinSearch(alpha, epsilon, g_vec)
    mean_arg = g_vec - alpha;
    mean_arg(mean_arg <= 0) = 0;
    y = 1 / epsilon * mean(mean_arg) + alpha;
end
