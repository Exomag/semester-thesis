function params = loadParameters()
    % Problem parameters
    params.HORIZON_TIME = 2;
    params.SAMPLING_TIME = 0.05;
    params.HORIZON = params.HORIZON_TIME / params.SAMPLING_TIME;
    params.CAR_WIDTH = 1;
    params.CAR_LENGTH = 2;
    params.LANE_WIDTH = 1.5;
    params.EGO_X0 = 0;
    params.EGO_Y0 = 1 / 2 * params.LANE_WIDTH;
    params.EGO_XDOT0 = 20;
    params.EGO_YDOT0 = 0;
    params.MIN_SPEED = 0;
    params.MAX_SPEED = 40;
    params.MIN_ACCEL = -40;
    params.MAX_ACCEL = 40;
    params.OBSTACLE_X = 10;
    params.OBSTACLE_Y = params.LANE_WIDTH / 2;
    params.ADV_X0 = 40;
    params.ADV_Y0 = 3 / 2 * params.LANE_WIDTH;
    params.ADV_XDOT0 = -20;
    params.ADV_YDOT0 = 0;

    % Dynamics
    params.A = [1, params.SAMPLING_TIME, 0, 0; 0, 1, 0, 0; 0, 0, 1, params.SAMPLING_TIME; 0, 0, 0, 1];
    params.B = [0, 0; params.SAMPLING_TIME, 0; 0, 0; 0, params.SAMPLING_TIME];

    % Objective parameters
    params.Q = diag([0, 10, 100, 10]);
    params.R = diag([1, 1]);
    params.M = 2000;

    % Chance constraint parameters
    params.BETA = 10^-3;
    params.EPSILON = 0.01;
    params.DIST_MEAN = params.ADV_XDOT0;
    params.DIST_STD = 5;
    params.DIST_VAR = params.DIST_STD^2;

    % Number of samples
    params.N_s = ceil(exp(1)/(exp(1) - 1)/params.EPSILON*(log((2^2) / params.BETA) + 2 * params.HORIZON - 1));

    % Risk allocation
    params.eps_single = params.EPSILON / params.HORIZON;
end
