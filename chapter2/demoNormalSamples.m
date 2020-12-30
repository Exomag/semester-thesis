function demoNormalSamples()
    % Parameters
    rng(0)
    beta = 10^-3;
    epsilon = 0.1;
    N_vals = (50:10:10000)';

    % Distribution parameters
    dist_mean = -0.5;
    dist_var = 1 / 12;
    dist_type = 'norm';
    x_star = norminv(1-epsilon) * sqrt(dist_var) + dist_mean;

    % Samples
    rng(0)
    samples = generateSamples(N_vals(end), dist_mean, dist_var, dist_type);
    N_scen = ceil(2/epsilon*log(1 / beta));

    % Trials
    x_scen = zeros(length(N_vals), 1);
    x_dist_known = zeros(length(N_vals), 1);
    x_cvar = zeros(length(N_vals), 1);
    for i = 1:length(N_vals)
        N = N_vals(i);

        % Scenario approach
        if N >= N_scen
            x_scen(i) = ScenarioApproach(samples(1:N));
        else
            x_scen(i) = Inf;
        end

        % Moments Robust approach
        x_dist_known(i) = KnownDist(samples(1:N), epsilon, beta);

        % CVaR
        x_cvar(i) = CvarApproach(samples(1:N), epsilon);
    end

    % Comparison plot
    figure();
    hold on;
    grid
    plot(N_vals, x_star*ones(size(N_vals)), 'LineWidth', 1)
    plot(N_vals, x_scen, 'LineWidth', 1)
    plot(N_vals, x_dist_known, 'LineWidth', 1)
    plot(N_vals, x_cvar, 'LineWidth', 1)
    xlabel('number of samples $$N$$', 'interpreter', 'latex')
    ylabel('optimization result', 'interpreter', 'latex')
    legend({'$$x^*$$', 'SA', 'MRA', 'CVaRA'}, 'interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    save2tikz('plots/CCProg_Norm_SamplesComp')
end

function x = ScenarioApproach(samples)
    f = 1;
    A = -ones(size(samples));
    b = -samples;
    options = optimoptions('linprog', 'Display', 'off');
    x = linprog(f, A, b, [], [], [], [], options);
end

function x = KnownDist(samples, epsilon, beta)
    beta = beta / 2;
    N = length(samples);
    mean_est = mean(samples);
    std_est = std(samples);
    var_est = var(samples);
    r1 = tinv(1-beta/2, N-1) * std_est / sqrt(N);
    r2 = var_est * max(abs(((N - 1)/chi2inv(beta / 2, N - 1)-1)), abs(abs(((N - 1)/chi2inv(1 - beta / 2, N - 1)-1))));
    fun = @(x) x;
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e3, 'Display', 'off');
    x = fmincon(fun, 0, [], [], [], [], [], [], @(x)KnownDistConstr(x, epsilon, mean_est, var_est, r1, r2), options);
end

function x = CvarApproach(samples, epsilon)
    fun = @(x) x;
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e3, 'Display', 'off');
    x = fmincon(fun, 0, [], [], [], [], [], [], @(x)CVarConstr(x, epsilon, samples), options);
end

function [c, ceq] = KnownDistConstr(x, epsilon, mean_est, var_est, r1, r2)
    c = norminv(1-epsilon) * sqrt(var_est+r2) - x + mean_est + r1;
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
