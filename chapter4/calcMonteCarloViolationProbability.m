function Pviol_joint = calcMonteCarloViolationProbability(params, solData)
    % Parameters
    HORIZON = params.HORIZON;
    CAR_WIDTH = params.CAR_WIDTH;
    CAR_LENGTH = params.CAR_LENGTH;
    ADV_X0 = params.ADV_X0;
    ADV_Y0 = params.ADV_Y0;
    ADV_YDOT0 = params.ADV_YDOT0;
    A = params.A;
    DIST_MEAN = params.DIST_MEAN;
    DIST_STD = params.DIST_STD;

    % Solution
    x_val = solData.x_val;

    % Initialization
    Nviol_joint = 0;
    rng(0)
    x_MonteCarlo = zeros(4, HORIZON+1);

    % Calculate violation probability via Monte Carlo sampling
    N_MonteCarlo = 1e5;
    tol = 1e-6;
    for i = 1:N_MonteCarlo
        sample_MonteCarlo = DIST_MEAN + DIST_STD * randn;
        x_MonteCarlo(:, 1) = [ADV_X0; sample_MonteCarlo; ADV_Y0; ADV_YDOT0];
        for j = 2:HORIZON + 1
            x_MonteCarlo(:, j) = A * x_MonteCarlo(:, j-1);
            if ~((x_val(1, j) + CAR_LENGTH / 2 - tol <= x_MonteCarlo(1, j) - CAR_LENGTH / 2) ...
                    || (x_val(1, j) - CAR_LENGTH / 2 + tol >= x_MonteCarlo(1, j) + CAR_LENGTH / 2) ...
                    || (x_val(3, j) + CAR_WIDTH / 2 - tol <= x_MonteCarlo(3, j) - CAR_WIDTH / 2) ...
                    || (x_val(3, j) - CAR_WIDTH / 2 + tol >= x_MonteCarlo(3, j) + CAR_WIDTH / 2))
                Nviol_joint = Nviol_joint + 1;
                break;
            end
        end
    end
    Pviol_joint = Nviol_joint / N_MonteCarlo;
end
