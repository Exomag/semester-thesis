function params = loadParameters()
    % Simulation
    params.HORIZON_TIME = 2;
    params.SAMPLING_TIME = 0.05;
    params.HORIZON = params.HORIZON_TIME / params.SAMPLING_TIME;

    % Dimensions
    params.CAR_WIDTH = 1;
    params.CAR_LENGTH = 2;
    params.LANE_WIDTH = 1.5;

    % Speed/Acceleration limit
    params.MIN_SPEED = 0;
    params.MAX_SPEED = 40;
    params.MIN_ACCEL = [-40; -40];
    params.MAX_ACCEL = [40; 40];

    % Ego initial state
    params.EGO_X0 = 0;
    params.EGO_Y0 = 1 / 2 * params.LANE_WIDTH;
    params.EGO_XDOT0 = 20;
    params.EGO_YDOT0 = 0;

    % Adversary 1 initial state
    params.ADV1_X0 = 10;
    params.ADV1_Y0 = 1 / 2 * params.LANE_WIDTH;
    params.ADV1_XDOT0 = 20;
    params.ADV1_YDOT0 = 0;

    % Adversary 2 initial state
    params.ADV2_X0 = 60;
    params.ADV2_Y0 = 3 / 2 * params.LANE_WIDTH;
    params.ADV2_XDOT0 = -20;
    params.ADV2_YDOT0 = 0;

    % Dynamics
    params.A = [1, params.SAMPLING_TIME, 0, 0; 0, 1, 0, 0; 0, 0, 1, params.SAMPLING_TIME; 0, 0, 0, 1];
    params.B = [0, 0; params.SAMPLING_TIME, 0; 0, 0; 0, params.SAMPLING_TIME];
    params.A_1D = params.A(1:2, 1:2);

    % Objective parameters
    params.Q = diag([0, 10, 5000, 0]);
    params.R = diag([1, 1]);
    params.P = 10000;
    params.M = 5000;

    % Chance constraint parameters
    params.BETA = 10^-3;
    params.EPSILON = 0.01;

    % Adversary 1 - x position
    params.DIST_ADV1_MEAN = 0;
    params.DIST_ADV1_XSTD = 0.1;
    params.DIST_ADV1_XVAR = params.DIST_ADV1_XSTD^2;

    % Adversary 1 - y position
    params.DIST_ADV1_YSTD = 0.05;
    params.DIST_ADV1_YVAR = params.DIST_ADV1_YSTD^2;

    % Adversary 2 - initial velocity
    params.DIST_ADV2_INIT_MEAN = [params.ADV2_X0; params.ADV2_XDOT0];
    params.DIST_ADV2_INIT_STD = diag([1, 10]);
    params.DIST_ADV2_INIT_VAR = params.DIST_ADV2_INIT_STD.^2;

    % Adversary 2 - process noise
    params.DIST_ADV2_NOISE_MEAN = [0; 0];
    params.DIST_ADV2_NOISE_STD = diag([0, 0.1]);
    params.DIST_ADV2_NOISE_VAR = params.DIST_ADV2_NOISE_STD.^2;

    % Number of samples
    params.NUMBER_SAMPLES = 1000;

    % Receding horizon parameters
    params.stepsHorizon = params.HORIZON;

    % Kalman Filter
    params.H = [1, 0];
    params.DIST_MEAS_NOISE_MEAN = 0;
    params.DIST_MEAS_NOISE_STD = 1;
    params.DIST_MEAS_NOISE_VAR = params.DIST_MEAS_NOISE_STD^2;
end
