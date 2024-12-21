% VAR Model Simulation and Estimation
clear;
clc;

% Parameters of the VAR model
T_values = [50, 100, 500, 1000, 5000]; % Different sample sizes
phi_true = [0.5, -0.2; 1.0, 0.8; -0.3, 0.5]; % True coefficients: [c1 c2; phi11 phi21; phi12 phi22]
sigma_u = [1.0, 0.5; 0.5, 1.0];             % Covariance matrix of errors

% Number of variables and lags
k = 2; % Number of variables in the VAR
p = 3; % Number of regressors (including intercept)

% Store results
phi_estimates = zeros(length(T_values), k * p);
phi_errors = zeros(length(T_values), k * p);

% Loop over different sample sizes
for t_idx = 1:length(T_values)
    T = T_values(t_idx);
    
    % Generate lagged dependent variables as regressors
    Y = zeros(T, k); % Placeholder for dependent variables
    X = zeros(T, p); % Placeholder for lagged variables and intercept

    % Initialize with random data
    Y(1:p, :) = randn(p, k); % Initial values for lagged variables

    for t = p+1:T
        % Construct the lagged matrix for current time
        X_t = [1, Y(t-1, 1), Y(t-1, 2)]; % Intercept + lagged values
        
        % Generate new values for Y
        Y(t, :) = X_t * phi_true + mvnrnd([0, 0], sigma_u);
    end

    % Build the full regressor matrix
    X = [ones(T, 1), Y(1:end-1, :)]; % Include intercept and lagged values

    % Stack the system
    Y_vec = reshape(Y(p+1:end, :), (T-p) * k, 1); % Stacked dependent variable (T*k x 1)
    X_block = kron(eye(k), X(p+1:end, :));        % Stacked regressor matrix (T*k x k*p)

    % Ordinary Least Squares (OLS) for initialization
    phi_ols = (X_block' * X_block) \ (X_block' * Y_vec);

    % Residual calculation
    residuals = Y_vec - X_block * phi_ols;
    residual_matrix = reshape(residuals, T-p, k); % Reshape residuals into (T-p) x k

    % Estimate residual covariance matrix
    Sigma_u_hat = (residual_matrix' * residual_matrix) / (T - p);

    % Generalized Least Squares (GLS) Estimation
    Sigma_u_inv = inv(Sigma_u_hat);
    W = kron(Sigma_u_inv, eye(T - p)); % Weighting matrix for GLS
    phi_gls = (X_block' * W * X_block) \ (X_block' * W * Y_vec);

    % Save results
    phi_estimates(t_idx, :) = phi_gls';
    phi_errors(t_idx, :) = phi_gls' - phi_true(:)';
end

% Display results
disp('True Parameters:');
disp(phi_true);

disp('Estimated Parameters and Errors for Different T:');
for t_idx = 1:length(T_values)
    T = T_values(t_idx);
    fprintf('T = %d\n', T);
    fprintf('Estimated Parameters:\n');
    disp(reshape(phi_estimates(t_idx, :), p, k));
    fprintf('Errors:\n');
    disp(reshape(phi_errors(t_idx, :), p, k));
end

% Plot convergence
figure;
for i = 1:numel(phi_true)
    plot(T_values, phi_errors(:, i), '-o', 'LineWidth', 1.5);
    hold on;
end
xlabel('Sample Size (T)');
ylabel('Parameter Estimation Error');
title('Convergence of Parameter Estimates');
legend('Phi11', 'Phi12', 'Phi13', 'Phi21', 'Phi22', 'Phi23');
grid on;
