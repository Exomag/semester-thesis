function [solData, solYALMIP, objective] = solveProblemSA(params)
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
    OBSTACLE_X = params.OBSTACLE_X;
    OBSTACLE_Y = params.OBSTACLE_Y;
    A = params.A;
    B = params.B;
    Q = params.Q;
    R = params.R;
    M = params.M;

    % Optimization parameters
    paramOpt = calcOptimizationParameters(params);
    x_adv_nom = paramOpt.x_adv_nom;
    realizations_advX = paramOpt.realizations_advX;
    x_0 = paramOpt.x_0;
    x_r = paramOpt.x_r;
    u_r = paramOpt.u_r;

    % Variables
    x_var = sdpvar(4, HORIZON+1);
    u_var = sdpvar(2, HORIZON);
    t_obs = binvar(4, HORIZON);
    t_adv = binvar(4, HORIZON);

    % Objective
    objective = 0;
    for i = 1:HORIZON
        objective = objective + (x_var(:, i + 1) - x_r)' * Q * (x_var(:, i + 1) - x_r) + (u_var(:, i) - u_r)' * R * (u_var(:, i) - u_r);
    end

    % Equality constraints
    constraints = (x_var(:, 1) == x_0);
    for i = 1:HORIZON
        constraints = [constraints, x_var(:, i + 1) == A * x_var(:, i) + B * u_var(:, i)];
    end

    % Inequality constraints
    for i = 1:HORIZON
        % Lanes
        constraints = [constraints, x_var(3, i + 1) - CAR_WIDTH / 2 >= 0];
        constraints = [constraints, x_var(3, i + 1) + CAR_WIDTH / 2 <= 2 * LANE_WIDTH];

        % Speed
        constraints = [constraints, MIN_SPEED <= x_var(2, i + 1) <= MAX_SPEED];

        % Acceleration
        constraints = [constraints, MIN_ACCEL <= u_var(:, i) <= MAX_ACCEL];

        % Obstacle
        constraints = [constraints, x_var(1, i + 1) + CAR_LENGTH / 2 <= OBSTACLE_X - CAR_LENGTH / 2 + M * t_obs(1, i)];
        constraints = [constraints, x_var(1, i + 1) - CAR_LENGTH / 2 >= OBSTACLE_X + CAR_LENGTH / 2 - M * t_obs(2, i)];
        constraints = [constraints, x_var(3, i + 1) + CAR_WIDTH / 2 <= OBSTACLE_Y - CAR_WIDTH / 2 + M * t_obs(3, i)];
        constraints = [constraints, x_var(3, i + 1) - CAR_WIDTH / 2 >= OBSTACLE_Y + CAR_WIDTH / 2 - M * t_obs(4, i)];
        constraints = [constraints, sum(t_obs(:, i)) <= 3];

        % Adversary
        constraints = [constraints, x_var(1, i + 1) + CAR_LENGTH / 2 <= realizations_advX(:, i + 1) - CAR_LENGTH / 2 + M * t_adv(1, i)];
        constraints = [constraints, x_var(1, i + 1) - CAR_LENGTH / 2 >= realizations_advX(:, i + 1) + CAR_LENGTH / 2 - M * t_adv(2, i)];
        constraints = [constraints, x_var(3, i + 1) + CAR_WIDTH / 2 <= x_adv_nom(3, i + 1) - CAR_WIDTH / 2 + M * t_adv(3, i)];
        constraints = [constraints, x_var(3, i + 1) - CAR_WIDTH / 2 >= x_adv_nom(3, i + 1) + CAR_WIDTH / 2 - M * t_adv(4, i)];
        constraints = [constraints, sum(t_adv(:, i)) <= 3];
    end

    % Solve
    options = sdpsettings('verbose', 0);
    solYALMIP = optimize(constraints, objective, options);

    % Solution data
    if solYALMIP.problem == 0
        solData.t_axis = SAMPLING_TIME * (0:HORIZON);
        solData.x_val = value(x_var);
        solData.u_val = value(u_var);
        solData.t_obs_val = value(t_obs);
        solData.t_adv_val = value(t_adv);
    else
        solData = [];
        disp('Hmm, something went wrong!');
        yalmiperror(solYALMIP.problem)
        return
    end
end
