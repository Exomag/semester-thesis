function [solData, solYALMIP] = solveProblemMRA(params, paramOpt)
    % Problem parameters
    SAMPLING_TIME = params.SAMPLING_TIME;
    HORIZON = params.HORIZON;
    CAR_WIDTH = params.CAR_WIDTH;
    CAR_LENGTH = params.CAR_LENGTH;
    LANE_WIDTH = params.LANE_WIDTH;
    MIN_SPEED = params.MIN_SPEED;
    MAX_SPEED = params.MAX_SPEED;
    MIN_ACCEL = params.MIN_ACCEL;
    MAX_ACCEL = params.MAX_ACCEL;
    A = params.A;
    B = params.B;
    Q = params.Q;
    R = params.R;
    P = params.P;
    M = params.M;

    % Optimization parameters
    adv2_NomState = paramOpt.adv2_NomState;
    adv1_Xmean = paramOpt.adv1_Xmean;
    adv1_Xvar = paramOpt.adv1_Xvar;
    adv1_Ymean = paramOpt.adv1_Ymean;
    adv1_Yvar = paramOpt.adv1_Yvar;
    adv2_Xmean = paramOpt.adv2_Xmean;
    adv2_Xvar = paramOpt.adv2_Xvar;
    adv1_X_r1 = paramOpt.adv1_X_r1;
    adv1_X_r2 = paramOpt.adv1_X_r2;
    adv1_Y_r1 = paramOpt.adv1_Y_r1;
    x_0 = paramOpt.x_0;
    x_r = paramOpt.x_r;
    u_r = paramOpt.u_r;
    risk = paramOpt.risk;

    %% Variables
    x_var = sdpvar(4, HORIZON+1);
    u_var = sdpvar(2, HORIZON);
    t_adv1 = binvar(4, HORIZON);
    t_adv2 = binvar(4, HORIZON);

    %% Objective
    objective = 0;
    for i = 1:HORIZON
        objective = objective + (x_var(:, i + 1) - x_r(:, i + 1))' * Q * (x_var(:, i + 1) - x_r(:, i + 1)) + (u_var(:, i) - u_r)' * R * (u_var(:, i) - u_r);
        objective = objective + max(P*(adv1_Xmean(i + 1) - x_var(1, i + 1)), 0);
    end

    %% Equality constraints
    constraints = x_var(:, 1) == x_0;
    for i = 1:HORIZON
        constraints = [constraints, x_var(:, i + 1) == A * x_var(:, i) + B * u_var(:, i)];
    end

    %% Inequality constraints
    for i = 1:HORIZON
        % Lanes
        constraints = [constraints, x_var(3, i + 1) - CAR_WIDTH / 2 >= 0];
        constraints = [constraints, x_var(3, i + 1) + CAR_WIDTH / 2 <= 2 * LANE_WIDTH];

        % Speed
        constraints = [constraints, MIN_SPEED <= x_var(2, i + 1) <= MAX_SPEED];

        % Acceleration
        constraints = [constraints, MIN_ACCEL <= u_var(:, i) <= MAX_ACCEL];

        % Adversary 1
        adv1_Xleq = -norminv(1-risk(1, i)) * sqrt(adv1_Xvar(i + 1)+adv1_X_r2) + adv1_Xmean(i+1) - adv1_X_r1;
        adv1_Xgeq = norminv(1-risk(2, i)) * sqrt(adv1_Xvar(i + 1)) + adv1_Xmean(i+1) + adv1_X_r1;
        adv1_Yleq = -norminv(1-risk(3, i)) * sqrt(adv1_Yvar(i + 1)) + adv1_Ymean(i+1) - adv1_Y_r1;
        adv1_Ygeq = norminv(1-risk(4, i)) * sqrt(adv1_Yvar(i + 1)) + adv1_Ymean(i+1) + adv1_Y_r1;
        constraints = [constraints, x_var(1, i + 1) + CAR_LENGTH / 2 <= adv1_Xleq - CAR_LENGTH / 2 + M * t_adv1(1, i)];
        constraints = [constraints, x_var(1, i + 1) - CAR_LENGTH / 2 >= adv1_Xgeq + CAR_LENGTH / 2 - M * t_adv1(2, i)];
        constraints = [constraints, x_var(3, i + 1) + CAR_WIDTH / 2 <= adv1_Yleq - CAR_WIDTH / 2 + M * t_adv1(3, i)];
        constraints = [constraints, x_var(3, i + 1) - CAR_WIDTH / 2 >= adv1_Ygeq + CAR_WIDTH / 2 - M * t_adv1(4, i)];
        constraints = [constraints, sum(t_adv1(:, i)) <= 3];

        % Adversary 2
        advX_leq = -norminv(1-risk(6, i)) * sqrt(adv2_Xvar(i + 1)) + adv2_Xmean(i+1);
        advX_geq = norminv(1-risk(5, i)) * sqrt(adv2_Xvar(i + 1)) + adv2_Xmean(i+1);
        constraints = [constraints, x_var(1, i + 1) + CAR_LENGTH / 2 <= advX_leq - CAR_LENGTH / 2 + M * t_adv2(1, i)];
        constraints = [constraints, x_var(1, i + 1) - CAR_LENGTH / 2 >= advX_geq + CAR_LENGTH / 2 - M * t_adv2(2, i)];
        constraints = [constraints, x_var(3, i + 1) + CAR_WIDTH / 2 <= adv2_NomState(3, i + 1) - CAR_WIDTH / 2 + M * t_adv2(3, i)];
        constraints = [constraints, x_var(3, i + 1) - CAR_WIDTH / 2 >= adv2_NomState(3, i + 1) + CAR_WIDTH / 2 - M * t_adv2(4, i)];
        constraints = [constraints, sum(t_adv2(:, i)) <= 3];
    end

    %% Solve
    options = sdpsettings('verbose', 0);
    solYALMIP = optimize(constraints, objective, options);

    %% Solution data
    if solYALMIP.problem == 0
        solData.t_axis = SAMPLING_TIME * (0:HORIZON);
        solData.x_val = value(x_var);
        solData.u_val = value(u_var);
        solData.t_adv1_val = value(t_adv1);
        solData.t_adv2_val = value(t_adv2);
    else
        solData = [];
    end
end
