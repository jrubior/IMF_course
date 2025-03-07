close all
clear
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

draws=cell([n_draws,3]);

for i=1:n_draws

    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*n,1) + reshape(PpsiTilde,n*m,1);
    Bdraw         = reshape(Bdraw,n*p+nex,n);
    draws{i,1}    = Bdraw;
    draws{i,2}    = Sigmadraw;
    draws{i,3}    = DrawQ(n);

end

function [Q,R] = DrawQ(n)
X=randn(n,n);
[Q,R]=qr(X);
for i=1:n
    if R(i,i) <0
        Q(:,i)=-Q(:,i);
        R(i,:)=-R(i,:);
    end
end
end

