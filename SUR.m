% SUR Model Simulation and Estimation
clear;
clc;

% Parameters of the SUR model
T_values = [50, 100, 500, 1000, 5000]; % Different sample sizes
beta_true = [0.5, 1.0; -0.2, 0.8];    % True coefficients: [beta11 beta12; beta21 beta22]
sigma_u = [1.0, 0.5; 0.5, 1.0];       % Covariance matrix of errors

% Number of equations and regressors
k = 2; % Number of equations
p = 2; % Number of regressors (including intercept)

% Store results
beta_estimates = zeros(length(T_values), k * p);
beta_errors = zeros(length(T_values), k * p);

% Loop over different sample sizes
for t_idx = 1:length(T_values)
    T = T_values(t_idx);
    
    % Generate regressors
    X = [ones(T, 1), randn(T, 1)]; % First column is intercept, second column is random

    % Generate error terms with given covariance matrix
    U = mvnrnd([0, 0], sigma_u, T); % T x 2 matrix of residuals
    
    % Generate dependent variables
    Y = zeros(T, k);
    for i = 1:T
        Y(i, 1) = X(i, :) * beta_true(:, 1) + U(i, 1);
        Y(i, 2) = X(i, :) * beta_true(:, 2) + U(i, 2);
    end
    
    % OLS estimation for each equation
    beta_ols = zeros(p, k);
    residuals = zeros(T, k);
    for eq = 1:k
        beta_ols(:, eq) = (X' * X) \ (X' * Y(:, eq));
        residuals(:, eq) = Y(:, eq) - X * beta_ols(:, eq);
    end
    
    % Estimate residual covariance matrix
    Sigma_u_hat = (residuals' * residuals) / T;

    % Generalized Least Squares (GLS) Estimation
    % Reshape system: Y (T*k x 1), X (T*k x k*p)
    Y_vec = reshape(Y, T * k, 1);
    X_block = kron(eye(k), X);
    Sigma_u_inv = inv(Sigma_u_hat);
    W = kron(Sigma_u_inv, eye(T)); % Weighting matrix
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
legend('Beta11', 'Beta12', 'Beta21', 'Beta22');
grid on;
