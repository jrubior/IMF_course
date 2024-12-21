% SUR Model Simulation and Estimation
clear;
clc;

% Parameters of the SUR model
T_values = [50, 100, 500, 1000, 5000]; % Different sample sizes
beta_true = [0.5, -0.2; 1.0, 0.8; -0.3, 0.5]; % True coefficients: [c1 c2; beta11 beta21; beta12 beta22]
sigma_u = [1.0, 0.5; 0.5, 1.0];       % Covariance matrix of errors

% Number of equations and regressors
k = 2; % Number of equations
p = 3; % Number of regressors (including intercept)

% Store results
beta_estimates = zeros(length(T_values), k * p);
beta_errors = zeros(length(T_values), k * p);

% Loop over different sample sizes
for t_idx = 1:length(T_values)
    T = T_values(t_idx);
    
    % Generate regressors
    X = [ones(T, 1), randn(T, 1), randn(T, 1)]; % Intercept + 2 random regressors

    % Generate error terms with given covariance matrix
    U = mvnrnd([0, 0], sigma_u, T); % T x 2 matrix of residuals
    
    % Generate dependent variables
    Y = zeros(T, k);
    for eq = 1:k
        Y(:, eq) = X * beta_true(:, eq) + U(:, eq); % Generate data for each equation
    end

    % Stack the system
    Y_vec = reshape(Y, T * k, 1); % Stacked dependent variable (T*k x 1)
    X_block = kron(eye(k), X);   % Stacked regressor matrix (T*k x k*p)

    % Ordinary Least Squares (OLS) for initialization
    beta_ols = (X_block' * X_block) \ (X_block' * Y_vec);

    % Residual calculation
    residuals = Y_vec - X_block * beta_ols;
    residual_matrix = reshape(residuals, T, k); % Reshape residuals into T x k

    % Estimate residual covariance matrix
    Sigma_u_hat = (residual_matrix' * residual_matrix) / T;

    % Generalized Least Squares (GLS) Estimation
    Sigma_u_inv = inv(Sigma_u_hat);
    W = kron(Sigma_u_inv, eye(T)); % Weighting matrix for GLS
    beta_gls = (X_block' * W * X_block) \ (X_block' * W * Y_vec);

    % Save results
    beta_estimates(t_idx, :) = beta_gls';
    beta_errors(t_idx, :) = beta_gls' - beta_true(:)';
end

% Display results
disp('True Parameters:');
disp(beta_true);

disp('Estimated Parameters and Errors for Different T:');
for t_idx = 1:length(T_values)
    T = T_values(t_idx);
    fprintf('T = %d\n', T);
    fprintf('Estimated Parameters:\n');
    disp(reshape(beta_estimates(t_idx, :), p, k));
    fprintf('Errors:\n');
    disp(reshape(beta_errors(t_idx, :), p, k));
end

% Plot convergence
figure;
for i = 1:numel(beta_true)
    plot(T_values, beta_errors(:, i), '-o', 'LineWidth', 1.5);
    hold on;
end
xlabel('Sample Size (T)');
ylabel('Parameter Estimation Error');
title('Convergence of Parameter Estimates');
legend('Beta11', 'Beta12', 'Beta13', 'Beta21', 'Beta22', 'Beta23');
grid on;

