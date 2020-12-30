function params = loadParameters()
    % Problem
    params.Ts = 0.5;
    params.c_x = 1;
    params.c_y = 2;
    params.l_x = 2;
    params.dv = 1;
    params.x1_0 = 0;
    params.x2_0 = 0;
    params.y1_0 = -params.l_x;
    params.y2_0 = -(params.c_y + 2);
    params.u_min = -2;
    params.u_max = 2;
    params.N = (4 + 2 * params.c_y) / params.Ts;

    % Disturbance
    params.beta = 10^-3;
    params.epsilon = 0.01;
    params.N_s = ceil(exp(1)/(exp(1) - 1)/params.epsilon*(log((2^(params.N)) / params.beta) + params.N * 1 - 1));
    params.dist_type = 'norm';
    params.dist_mean = 0;
    params.dist_std = abs(params.y1_0) / 6;
    params.dist_var = params.dist_std^2;

    % Dynamics
    params.A = [1, params.Ts; 0, 1];
    params.B = [0; params.Ts];

    % Optimization data
    params.nx = 2;
    params.nu = 1;
    params.Q = diag([10, 1]);
    params.R = 1;
    params.M = 100;
end
