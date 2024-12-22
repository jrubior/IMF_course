close all
clear all
clc
T = 1000;
n_draws=100000;

beta_true = [0.1, 0.2; 0.9, 0.08; -0.03, 0.75]; % True coefficients: [c1 c2; beta11 beta21; beta12 beta22]
sigma_u = [1.0, 0.5; 0.5, 1.0];       % Covariance matrix of errors


p   = 1;            % number of lags
n   = 2;            % number of endogenous variables
nex = 1;            % set equal to 1 if a constant is included; 0 otherwise
m   = n*p + nex;    % number of independent variables
e  = eye(n);        % create identity matrix

Y = zeros(T, n);

Y(1,:) = [randn(1, 1), randn(1, 1)];

% Generate error terms with given covariance matrix

U = mvnrnd([0, 0], sigma_u, T); % T x 2 matrix of residuals

% Generate dependent variables

for t=2:T

    Y(t, :) = [1,Y(t-1,:)] * beta_true + U(t, :);

end

Ytemp=Y(2:end,:);
X=[Y(1:end-1,:) ones(T-1,1)]; % this is the notation used in the SUR
X=[ones(T-1,1) Y(1:end-1,:) ]; % this is the notation used in our papers
Y=Ytemp;


%% prior for reduced-form parameters
nnuBar              = 0;
OomegaBarInverse    = zeros(m);
PpsiBar             = zeros(m,n);
PphiBar             = zeros(n);

%% posterior for reduced-form parameters
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
PpsiTilde           = OomegaTilde*(X'*Y + OomegaBarInverse*PpsiBar);
PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - PpsiTilde'*OomegaTildeInverse*PpsiTilde;
PphiTilde           = (PphiTilde'+PphiTilde)*0.5;

hh       = @(x)chol(x);

cholOomegaTilde=hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

%% drawing Sigma and B|Sigma

draws=cell([n_draws,2]);

for i=1:n_draws

Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
cholSigmadraw = hh(Sigmadraw)';
Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*n,1) + reshape(PpsiTilde,n*m,1);
Bdraw         = reshape(Bdraw,n*p+nex,n);
draws{i,1}    = Bdraw;
draws{i,2}    = Sigmadraw;

end

around=.15;

figure;
for i = 1:2
    for j = 1:2
        subplot(2, 2, (i-1)*2 + j);
        sigma_draws = cellfun(@(x) x(i, j), draws(:, 2));
        histogram(sigma_draws, 30);
        hold on;
        xline(sigma_u(i, j), 'r', 'LineWidth', 2);
        xlim([sigma_u(i, j) - around, sigma_u(i, j) + around]);
        title(['\sigma_{' num2str(i) ',' num2str(j) '}']);
        hold off;
    end
end


figure;
i=1;
for j = 1:2
subplot(3, 2, (i-1)*2 + j);
beta_draws = cellfun(@(x) x(i, j), draws(:, 1));
histogram(beta_draws, 30);
hold on;
xline(beta_true(i, j), 'r', 'LineWidth', 2);
xlim([beta_true(i, j) - around, beta_true(i, j) + around]);
title(['c_{' num2str(j) '}']);
end
for i = 2:3
    for j = 1:2
        subplot(3, 2, (i-1)*2 + j);
        beta_draws = cellfun(@(x) x(i, j), draws(:, 1));
        histogram(beta_draws, 30);
        hold on;
        xline(beta_true(i, j), 'r', 'LineWidth', 2);
        xlim([beta_true(i, j) - around, beta_true(i, j) + around]);
        title(['\beta_{' num2str(i-1) ',' num2str(j) '}']);
        hold off;
    end
end