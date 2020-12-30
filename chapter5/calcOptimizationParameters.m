function paramOpt = calcOptimizationParameters(params)
    % Parameters
    HORIZON = params.HORIZON;
    EGO_X0 = params.EGO_X0;
    EGO_Y0 = params.EGO_Y0;
    EGO_XDOT0 = params.EGO_XDOT0;
    EGO_YDOT0 = params.EGO_YDOT0;
    ADV1_X0 = params.ADV1_X0;
    ADV1_Y0 = params.ADV1_Y0;
    ADV1_XDOT0 = params.ADV1_XDOT0;
    ADV1_YDOT0 = params.ADV1_YDOT0;
    ADV2_X0 = params.ADV2_X0;
    ADV2_Y0 = params.ADV2_Y0;
    ADV2_XDOT0 = params.ADV2_XDOT0;
    ADV2_YDOT0 = params.ADV2_YDOT0;
    A = params.A;
    A_1D = params.A_1D;
    BETA = params.BETA;
    EPSILON = params.EPSILON;
    DIST_ADV1_MEAN = params.DIST_ADV1_MEAN;
    DIST_ADV1_XSTD = params.DIST_ADV1_XSTD;
    DIST_ADV1_XVAR = params.DIST_ADV1_XVAR;
    DIST_ADV1_YSTD = params.DIST_ADV1_YSTD;
    DIST_ADV1_YVAR = params.DIST_ADV1_YVAR;
    DIST_ADV2_INIT_MEAN = params.DIST_ADV2_INIT_MEAN;
    DIST_ADV2_INIT_VAR = params.DIST_ADV2_INIT_VAR;
    DIST_ADV2_NOISE_MEAN = params.DIST_ADV2_NOISE_MEAN;
    DIST_ADV2_NOISE_VAR = params.DIST_ADV2_NOISE_VAR;
    NUMBER_SAMPLES = params.NUMBER_SAMPLES;

    % Adversary 1 nominal trajectory
    adv1_NomState = [ADV1_X0; ADV1_XDOT0; ADV1_Y0; ADV1_YDOT0];
    for i = 1:HORIZON
        adv1_NomState = [adv1_NomState, A * adv1_NomState(:, end)];
    end
    paramOpt.adv1_NomState = adv1_NomState;

    % Adversary 2 nominal trajectory
    adv2_NomState = [ADV2_X0; ADV2_XDOT0; ADV2_Y0; ADV2_YDOT0];
    for i = 1:HORIZON
        adv2_NomState = [adv2_NomState, A * adv2_NomState(:, end)];
    end
    paramOpt.adv2_NomState = adv2_NomState;

    % Adversary 1 moments
    paramOpt.adv1_Xmean = adv1_NomState(1, :);
    paramOpt.adv1_Xstd = DIST_ADV1_XSTD * ones(size(paramOpt.adv1_Xmean));
    paramOpt.adv1_Xvar = DIST_ADV1_XVAR * ones(size(paramOpt.adv1_Xmean));
    paramOpt.adv1_Ymean = adv1_NomState(3, :);
    paramOpt.adv1_Ystd = DIST_ADV1_YSTD * ones(size(paramOpt.adv1_Ymean));
    paramOpt.adv1_Yvar = DIST_ADV1_YVAR * ones(size(paramOpt.adv1_Ymean));

    % Adversary 2 moments
    adv2_Xmean = zeros(size(adv2_NomState(1, :)));
    adv2_Xstd = zeros(size(adv2_Xmean));
    adv2_Xvar = zeros(size(adv2_Xmean));
    for i = 1:length(adv2_Xmean)
        temp_mean = A_1D^(i - 1) * DIST_ADV2_INIT_MEAN;
        temp_var = A_1D^(i - 1) * DIST_ADV2_INIT_VAR * (A_1D^(i - 1))';
        for k = 0:i - 2
            temp_mean = temp_mean + A_1D^(i - k - 2) * DIST_ADV2_NOISE_MEAN;
            temp_var = temp_var + A_1D^k * DIST_ADV2_NOISE_VAR * (A_1D^k)';
        end
        adv2_Xmean(i) = temp_mean(1);
        adv2_Xstd(i) = sqrt(temp_var(1, 1));
        adv2_Xvar(i) = temp_var(1, 1);
    end
    paramOpt.adv2_Xmean = adv2_Xmean;
    paramOpt.adv2_Xstd = adv2_Xstd;
    paramOpt.adv2_Xvar = adv2_Xvar;

    % x position moment bounds
    paramOpt.adv1_Xsamples = DIST_ADV1_MEAN + DIST_ADV1_XSTD * randn(NUMBER_SAMPLES, 1);
    paramOpt.adv1_Xmean_est = mean(paramOpt.adv1_Xsamples);
    paramOpt.adv1_Xstd_est = std(paramOpt.adv1_Xsamples);
    paramOpt.adv1_Xvar_est = var(paramOpt.adv1_Xsamples);
    paramOpt.adv1_X_r1 = tinv(1-BETA/2, NUMBER_SAMPLES-1) .* paramOpt.adv1_Xstd_est ./ sqrt(NUMBER_SAMPLES);
    paramOpt.adv1_X_r2 = paramOpt.adv1_Xvar_est * max(abs(((NUMBER_SAMPLES - 1)/chi2inv(BETA / 2, NUMBER_SAMPLES - 1)-1)), abs(abs(((NUMBER_SAMPLES - 1)/chi2inv(1 - BETA / 2, NUMBER_SAMPLES - 1)-1))));

    % y position moment bounds
    paramOpt.adv1_Ysamples = DIST_ADV1_MEAN + DIST_ADV1_YSTD * randn(NUMBER_SAMPLES, 1);
    paramOpt.adv1_Ymean_est = mean(paramOpt.adv1_Ysamples);
    paramOpt.adv1_Ystd_est = std(paramOpt.adv1_Ysamples);
    paramOpt.adv1_Yvar_est = var(paramOpt.adv1_Ysamples);
    paramOpt.adv1_Y_r1 = tinv(1-BETA/2, NUMBER_SAMPLES-1) .* paramOpt.adv1_Ystd_est ./ sqrt(NUMBER_SAMPLES);
    paramOpt.adv1_Y_r2 = paramOpt.adv1_Yvar_est * max(abs(((NUMBER_SAMPLES - 1)/chi2inv(BETA / 2, NUMBER_SAMPLES - 1)-1)), abs(abs(((NUMBER_SAMPLES - 1)/chi2inv(1 - BETA / 2, NUMBER_SAMPLES - 1)-1))));

    % References
    paramOpt.x_0 = [EGO_X0; EGO_XDOT0; EGO_Y0; EGO_YDOT0];
    paramOpt.x_r = [EGO_X0 * ones(1, HORIZON + 1); EGO_XDOT0 * ones(1, HORIZON + 1); EGO_Y0 * ones(1, HORIZON + 1); EGO_YDOT0 * ones(1, HORIZON + 1)];
    paramOpt.u_r = [0; 0];

    % Risk allocation
    paramOpt.risk = EPSILON / (3 * HORIZON) * ones(6, HORIZON);
end
