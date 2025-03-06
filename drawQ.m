clear 
close all
clc

n=3;

Q = zeros(n, n);

for j = 1:n
    % Draw x_j independently from a standard normal distribution
    x_j = randn(n + 1 - j, 1);

    % Normalize to get w_j
    w_j = x_j / norm(x_j);

    % Compute M_j
    if j == 1
        K_j= eye(n);
    else
        M_j = Q(:, 1:(j-1))';
        K_j= null(M_j);
    end    

    % Compute q_j
    q_j = K_j * w_j;

    % Store q_j in Q
    Q(:, j) = q_j;
end
