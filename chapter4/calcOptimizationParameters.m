function paramOpt = calcOptimizationParameters(params)
    % Parameters
    HORIZON_TIME = params.HORIZON_TIME;
    HORIZON = params.HORIZON;
    EGO_X0 = params.EGO_X0;
    EGO_Y0 = params.EGO_Y0;
    EGO_XDOT0 = params.EGO_XDOT0;
    EGO_YDOT0 = params.EGO_YDOT0;
    ADV_X0 = params.ADV_X0;
    ADV_Y0 = params.ADV_Y0;
    ADV_XDOT0 = params.ADV_XDOT0;
    ADV_YDOT0 = params.ADV_YDOT0;
    A = params.A;
    BETA = params.BETA;
    DIST_MEAN = params.DIST_MEAN;
    DIST_STD = params.DIST_STD;
    N_s = params.N_s;

    % Adversary nominal trajectory
    x_adv_nom = [ADV_X0; ADV_XDOT0; ADV_Y0; ADV_YDOT0];
    for i = 1:HORIZON
        x_adv_nom = [x_adv_nom, A * x_adv_nom(:, end)];
    end
    paramOpt.x_adv_nom = x_adv_nom;

    % Disturbance realizations
    samples = DIST_MEAN + DIST_STD * randn(N_s, 1);
    realizations_advX = zeros(N_s, HORIZON+1);
    x_adv = zeros(4, HORIZON+1);
    for i = 1:N_s
        x_adv(:, 1) = [ADV_X0; samples(i); ADV_Y0; ADV_YDOT0];
        for j = 2:HORIZON + 1
            x_adv(:, j) = A * x_adv(:, j-1);
        end
        realizations_advX(i, :) = x_adv(1, :);
    end
    paramOpt.samples = samples;
    paramOpt.realizations_advX = realizations_advX;
    paramOpt.x_adv = x_adv;

    % Disturbance moments
    paramOpt.beta_mra = BETA / 2;
    paramOpt.mean_advX = mean(realizations_advX);
    paramOpt.var_advX = var(realizations_advX);
    paramOpt.std_advX = std(realizations_advX);
    paramOpt.r1 = tinv(1-paramOpt.beta_mra/2, N_s-1) .* paramOpt.std_advX ./ sqrt(N_s);
    paramOpt.r2 = paramOpt.var_advX * max(abs(((N_s - 1)/chi2inv(paramOpt.beta_mra / 2, N_s - 1)-1)), abs(abs(((N_s - 1)/chi2inv(1 - paramOpt.beta_mra / 2, N_s - 1)-1))));

    % References
    paramOpt.x_0 = [EGO_X0; EGO_XDOT0; EGO_Y0; EGO_YDOT0];
    paramOpt.x_r = [EGO_X0 + HORIZON_TIME * EGO_XDOT0; EGO_XDOT0; EGO_Y0; EGO_YDOT0];
    paramOpt.u_r = [0; 0];
end
