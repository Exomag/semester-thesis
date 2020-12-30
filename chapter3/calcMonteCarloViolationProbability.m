function [Pviol_joint, Pviol_single] = calcMonteCarloViolationProbability(params, solData)
    % Parameters
    c_x = params.c_x;
    N = params.N;
    dist_type = params.dist_type;
    dist_mean = params.dist_mean;
    dist_var = params.dist_var;
    M = params.M;

    % Solution
    t = solData.t;
    x1_abs = solData.x1_abs;
    y1_abs = solData.y1_abs;

    % Calculate violation probability via Monte Carlo sampling
    N_MonteCarlo = 1e6;
    Nviol_single = zeros(N, 1);
    Nviol_joint = 0;
    rng(0)
    for i = 1:N_MonteCarlo
        samples_trial = generateSamples(N, dist_mean, dist_var, dist_type);
        if any(x1_abs(2:end) < y1_abs(2:end)+c_x-M*t(:, 2)+samples_trial)
            Nviol_joint = Nviol_joint + 1;
        end
        Nviol_single = Nviol_single + double(x1_abs(2:end) < y1_abs(2:end)+c_x-M*t(:, 2)+samples_trial);
    end
    Pviol_joint = Nviol_joint / N_MonteCarlo;
    Pviol_single = Nviol_single / N_MonteCarlo;
end
