function paramOpt = updateOptimizationParameters(params, paramOpt, stateSim)
    % Parameters
    HORIZON = params.HORIZON;
    A = params.A;
    A_1D = params.A_1D;
    DIST_ADV1_XSTD = params.DIST_ADV1_XSTD;
    DIST_ADV1_XVAR = params.DIST_ADV1_XVAR;
    DIST_ADV1_YSTD = params.DIST_ADV1_YSTD;
    DIST_ADV1_YVAR = params.DIST_ADV1_YVAR;
    DIST_ADV2_NOISE_MEAN = params.DIST_ADV2_NOISE_MEAN;
    DIST_ADV2_NOISE_VAR = params.DIST_ADV2_NOISE_VAR;

    % Adversary 1 nominal trajectory
    adv1_NomState = stateSim.adv1_State;
    for i = 1:HORIZON
        adv1_NomState = [adv1_NomState, A * adv1_NomState(:, end)];
    end
    paramOpt.adv1_NomState = adv1_NomState;

    % Adversary 2 nominal trajectory
    adv2_NomState = stateSim.adv2_State;
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
        temp_mean = A_1D^(i - 1) * stateSim.xm;
        temp_var = A_1D^(i - 1) * stateSim.Pm * (A_1D^(i - 1))';
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

    % References
    paramOpt.x_0 = stateSim.ego_State;
end
