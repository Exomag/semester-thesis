function simData = solveProblemRecedingHorizon(params, demo_mode)
    % Problem parameters
    SAMPLING_TIME = params.SAMPLING_TIME;
    CAR_WIDTH = params.CAR_WIDTH;
    CAR_LENGTH = params.CAR_LENGTH;
    EGO_X0 = params.EGO_X0;
    EGO_Y0 = params.EGO_Y0;
    EGO_XDOT0 = params.EGO_XDOT0;
    EGO_YDOT0 = params.EGO_YDOT0;
    ADV1_X0 = params.ADV1_X0;
    ADV1_Y0 = params.ADV1_Y0;
    ADV1_XDOT0 = params.ADV1_XDOT0;
    ADV1_YDOT0 = params.ADV1_YDOT0;
    ADV2_Y0 = params.ADV2_Y0;
    ADV2_YDOT0 = params.ADV2_YDOT0;
    A = params.A;
    A_1D = params.A_1D;
    Q = params.Q;
    R = params.R;
    P = params.P;
    DIST_ADV1_MEAN = params.DIST_ADV1_MEAN;
    DIST_ADV1_XSTD = params.DIST_ADV1_XSTD;
    DIST_ADV1_YSTD = params.DIST_ADV1_YSTD;
    DIST_ADV2_INIT_MEAN = params.DIST_ADV2_INIT_MEAN;
    DIST_ADV2_INIT_VAR = params.DIST_ADV2_INIT_VAR;
    DIST_ADV2_NOISE_MEAN = params.DIST_ADV2_NOISE_MEAN;
    DIST_ADV2_NOISE_VAR = params.DIST_ADV2_NOISE_VAR;
    stepsHorizon = params.stepsHorizon;
    H = params.H;
    DIST_MEAS_NOISE_MEAN = params.DIST_MEAS_NOISE_MEAN;
    DIST_MEAS_NOISE_VAR = params.DIST_MEAS_NOISE_VAR;

    % Optimization parameters
    paramOpt = calcOptimizationParameters(params);
    adv1_Xmean = paramOpt.adv1_Xmean;
    x_r = paramOpt.x_r;
    u_r = paramOpt.u_r;
    % Optimization problem
    tHorizon = SAMPLING_TIME * (0:stepsHorizon);
    ego_StatePlan = cell(stepsHorizon+1, 1);
    simData.totalTime = 0;
    simData.totalCost = 0;
    simData.feasibleProblem = true;
    simData.collision = false;

    % States
    ego_State = zeros(4, stepsHorizon+1);
    ego_Input = zeros(2, stepsHorizon);
    adv1_State = zeros(4, stepsHorizon+1);
    adv1_StateNoisy = zeros(4, stepsHorizon+1);
    adv2_State = zeros(4, stepsHorizon+1);

    % Kalman filter
    xm = cell(stepsHorizon+1, 1);
    xp = cell(stepsHorizon, 1);
    K = cell(stepsHorizon, 1);
    Pm = cell(stepsHorizon+1, 1);
    Pp = cell(stepsHorizon, 1);
    adv2_Xmeas = zeros(stepsHorizon, 1);
    adv2_mean_est = zeros(2, stepsHorizon+1);
    adv2_var_est = zeros(2, stepsHorizon+1);

    % Initial states
    ego_State(:, 1) = [EGO_X0; EGO_XDOT0; EGO_Y0; EGO_YDOT0];
    adv1_State(:, 1) = [ADV1_X0; ADV1_XDOT0; ADV1_Y0; ADV1_YDOT0];
    adv1_StateNoisy(:, 1) = [ADV1_X0; ADV1_XDOT0; ADV1_Y0; ADV1_YDOT0];
    adv2_State(:, 1) = [mvnrnd(DIST_ADV2_INIT_MEAN, DIST_ADV2_INIT_VAR)'; ADV2_Y0; ADV2_YDOT0];
    if demo_mode
        adv2_State(:, 1) = [DIST_ADV2_INIT_MEAN + [2; 2]; ADV2_Y0; ADV2_YDOT0];
    end

    % Initial estimates
    xm{1} = DIST_ADV2_INIT_MEAN;
    Pm{1} = DIST_ADV2_INIT_VAR;
    adv2_mean_est(:, 1) = xm{1};
    adv2_var_est(:, 1) = [Pm{1}(1, 1); Pm{1}(2, 2)];

    % Receding Horizon
    for indHorizon = 1:stepsHorizon
        % Solve problem
        stateSim.adv1_State = adv1_State(:, indHorizon);
        stateSim.adv2_State = adv2_State(:, indHorizon);
        stateSim.xm = xm{indHorizon};
        stateSim.Pm = Pm{indHorizon};
        stateSim.ego_State = ego_State(:, indHorizon);
        paramOpt = updateOptimizationParameters(params, paramOpt, stateSim);
        [solData, solution] = solveProblemMRA(params, paramOpt);

        % Update time/cost
        if solution.problem == 0
            x_val = solData.x_val;
            u_val = solData.u_val;
            simData.totalTime = simData.totalTime + solution.solvertime;
            simData.totalCost = simData.totalCost + (x_val(:, 2) - x_r(:, 2))' * Q * (x_val(:, 2) - x_r(:, 2)) + (u_val(:, 1) - u_r)' * R * (u_val(:, 1) - u_r) + max(P*(adv1_Xmean(2) - x_val(1, 2)), 0);
        else
            yalmiperror(solution.problem)
            simData.feasibleProblem = false;
            return;
        end

        % Update states
        ego_State(:, indHorizon+1) = x_val(:, 2);
        ego_Input(:, indHorizon) = u_val(:, 1);
        adv1_State(:, indHorizon+1) = A * adv1_State(:, indHorizon);
        adv1_StateNoisy(:, indHorizon+1) = adv1_State(:, indHorizon+1) + [DIST_ADV1_MEAN + DIST_ADV1_XSTD * randn(); 0; DIST_ADV1_MEAN + DIST_ADV1_YSTD * randn(); 0];
        adv2_State(:, indHorizon+1) = A * adv2_State(:, indHorizon) + [mvnrnd(DIST_ADV2_NOISE_MEAN, DIST_ADV2_NOISE_VAR)'; zeros(2, 1)];
        adv2_Xmeas(indHorizon) = adv2_State(1, indHorizon+1) + mvnrnd(DIST_MEAS_NOISE_MEAN, DIST_MEAS_NOISE_VAR);
        ego_StatePlan{indHorizon} = x_val;

        % Kalman filter
        xp{indHorizon} = A_1D * xm{indHorizon};
        Pp{indHorizon} = A_1D * Pm{indHorizon} * A_1D' + DIST_ADV2_NOISE_VAR;
        K{indHorizon} = Pp{indHorizon} * H' * inv(H*Pp{indHorizon}*H'+DIST_MEAS_NOISE_VAR);
        xm{indHorizon+1} = xp{indHorizon} + K{indHorizon} * (adv2_Xmeas(indHorizon) - H * xp{indHorizon});
        Pm{indHorizon+1} = (eye(2) - K{indHorizon} * H) * Pp{indHorizon} * (eye(2) - K{indHorizon} * H)' + K{indHorizon} * DIST_MEAS_NOISE_VAR * K{indHorizon}';

        % Update estimates
        adv2_mean_est(:, indHorizon+1) = xm{indHorizon+1};
        adv2_var_est(:, indHorizon+1) = [Pm{indHorizon + 1}(1, 1); Pm{indHorizon + 1}(2, 2)];
    end

    % Finishing up
    indHorizon = indHorizon + 1;
    stateSim.adv1_State = adv1_State(:, indHorizon);
    stateSim.adv2_State = adv2_State(:, indHorizon);
    stateSim.xm = xm{indHorizon};
    stateSim.Pm = Pm{indHorizon};
    stateSim.ego_State = ego_State(:, indHorizon);
    paramOpt = updateOptimizationParameters(params, paramOpt, stateSim);
    [solData, solution] = solveProblemMRA(params, paramOpt);
    if solution.problem == 0
        x_val = solData.x_val;
    else
        yalmiperror(solution.problem)
        simData.feasibleProblem = false;
        return;
    end
    ego_StatePlan{indHorizon} = x_val;
    simData.adv1_StateNoNoise = adv1_State;
    adv1_State = adv1_StateNoisy;

    % Collision check
    tol = 1e-6;
    for i = 2:stepsHorizon + 1
        if ~((ego_State(1, i) + CAR_LENGTH / 2 - tol <= adv1_State(1, i) - CAR_LENGTH / 2) ...
                || (ego_State(1, i) - CAR_LENGTH / 2 + tol >= adv1_State(1, i) + CAR_LENGTH / 2) ...
                || (ego_State(3, i) + CAR_WIDTH / 2 - tol <= adv1_State(3, i) - CAR_WIDTH / 2) ...
                || (ego_State(3, i) - CAR_WIDTH / 2 + tol >= adv1_State(3, i) + CAR_WIDTH / 2)) ...
                || ~((ego_State(1, i) + CAR_LENGTH / 2 - tol <= adv2_State(1, i) - CAR_LENGTH / 2) ...
                || (ego_State(1, i) - CAR_LENGTH / 2 + tol >= adv2_State(1, i) + CAR_LENGTH / 2) ...
                || (ego_State(3, i) + CAR_WIDTH / 2 - tol <= adv2_State(3, i) - CAR_WIDTH / 2) ...
                || (ego_State(3, i) - CAR_WIDTH / 2 + tol >= adv2_State(3, i) + CAR_WIDTH / 2))
            simData.collision = true;
        end
    end

    % Output data
    simData.tHorizon = tHorizon;
    simData.ego_StatePlan = ego_StatePlan;
    simData.x_val = ego_State;
    simData.u_val = ego_Input;
    simData.adv1_NomState = adv1_State;
    simData.adv2_NomState = adv2_State;
    simData.adv2_mean_est = adv2_mean_est;
    simData.adv2_var_est = adv2_var_est;
end