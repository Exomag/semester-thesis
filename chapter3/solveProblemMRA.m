function [solData, solYALMIP, objective, samples, k_active] = solveProblemMRA(params)
    % Parameters
    Ts = params.Ts;
    c_x = params.c_x;
    c_y = params.c_y;
    l_x = params.l_x;
    dv = params.dv;
    x1_0 = params.x1_0;
    x2_0 = params.x2_0;
    y1_0 = params.y1_0;
    y2_0 = params.y2_0;
    u_min = params.u_min;
    u_max = params.u_max;
    N = params.N;
    beta = params.beta;
    epsilon = params.epsilon;
    N_s = params.N_s;
    dist_type = params.dist_type;
    dist_mean = params.dist_mean;
    dist_var = params.dist_var;
    A = params.A;
    B = params.B;
    nx = params.nx;
    nu = params.nu;
    Q = params.Q;
    R = params.R;
    M = params.M;

    % Samples
    beta_mra = beta / 2;
    for i = 1:N
        samples{i} = generateSamples(N_s, dist_mean, dist_var, dist_type);
        mean_est{i} = mean(samples{i});
        std_est{i} = std(samples{i});
        var_est{i} = var(samples{i});
        r1{i} = tinv(1-beta_mra/2, N_s-1) * std_est{i} / sqrt(N_s);
        r2{i} = var_est{i} * max(abs(((N_s - 1)/chi2inv(beta_mra / 2, N_s - 1)-1)), abs(abs(((N_s - 1)/chi2inv(1 - beta_mra / 2, N_s - 1)-1))));
    end

    % Uniform risk allocation
    eps_single = epsilon / N;

    % Car trajectory
    y1 = [y1_0; zeros(N, 1)];
    y2 = [y2_0; zeros(N, 1)];
    k_active = [];
    for i = 2:N + 1
        y1(i) = y1(i-1);
        y2(i) = y2(i-1) + dv * Ts;
        if (y2(i) + c_y / 2 >= -c_y / 2) && (y2(i) - c_y / 2 <= c_y / 2)
            k_active = [k_active; i - 1];
        end
    end

    % Variables
    u_var = sdpvar(nu*ones(1, N), ones(1, N));
    x_var = sdpvar(nx*ones(1, N + 1), ones(1, N + 1));
    t_var = binvar(4*ones(1, N), ones(1, N));

    % Constraints & Objective
    constraints = [];
    objective = 0;
    x0 = [x1_0; x2_0];
    constraints = [constraints, x_var{1} == x0];
    for k = 1:N
        x_var_curr = x_var{k+1};
        x1_var_curr = x_var_curr(1);
        objective = objective + norm(Q*x_var{k}, 2) + norm(R*u_var{k}, 2);
        constraints = [constraints, x_var{k + 1} == A * x_var{k} + B * u_var{k}];
        constraints = [constraints, u_min <= u_var{k} <= u_max, -l_x / 2 + c_x / 2 <= x1_var_curr <= l_x / 2 - c_x / 2];
        t_var_curr = t_var{k};
        constraints = [constraints, x1_var_curr <= y1(k + 1) - c_x + M * t_var_curr(1)];
        constraints = [constraints, x1_var_curr >= y1(k + 1) + c_x - M * t_var_curr(2) + (norminv(1 - eps_single) * sqrt(var_est{k} + r2{k}) + mean_est{k} + r1{k})];
        constraints = [constraints, M * t_var_curr(3) >= - y2(k + 1) + c_y];
        constraints = [constraints, M * t_var_curr(4) >= y2(k + 1) + c_y];
        constraints = [constraints, sum(t_var_curr) <= 3];
    end

    % Solve
    options = sdpsettings('verbose', 0);
    solYALMIP = optimize(constraints, objective, options);

    % Solution
    if solYALMIP.problem == 0
        x_curr = value(x_var{1});
        x1 = x_curr(1);
        x2 = x_curr(2);
        u = [];
        t = [];
        for k = 1:N
            x_curr = value(x_var{k + 1});
            x1 = [x1; x_curr(1)];
            x2 = [x2; x_curr(2)];
            u = [u; value(u_var{k})];
            t = [t; value(t_var{k})'];
        end
    else
        solData = [];
        disp('Hmm, something went wrong!');
        solYALMIP.info
        yalmiperror(solYALMIP.problem)
        return
    end

    % Absolute data
    v1 = 20;
    v2 = v1 + dv;
    x1_abs = x1 + l_x / 2;
    y1_abs = y1 + l_x / 2;
    x2_abs = [x2_0; zeros(N, 1)];
    y2_abs = [y2_0; zeros(N, 1)];
    for i = 2:N + 1
        x2_abs(i) = x2_abs(i-1) + v1 * Ts;
        y2_abs(i) = y2_abs(i-1) + v2 * Ts;
    end
    
    % Output data
    solData.t = t;
    solData.x1 = x1;
    solData.y1 = y1;
    solData.x2 = x2;
    solData.y2 = y2;
    solData.x1_abs = x1_abs;
    solData.y1_abs = y1_abs;
    solData.x2_abs = x2_abs;
    solData.y2_abs = y2_abs;
end
