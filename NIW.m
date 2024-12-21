T = 100;
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
X=[ones(T-1,1) Y(1:end-1,:)];
Y=Ytemp;

beta_ols = (X' * X) \ (X' * Y);
residuals= Y - X * beta_ols;


% Estimate residual covariance matrix
Sigma_u_hat = (residuals' * residuals) / T;

% Generalized Least Squares (GLS) Estimation
% Reshape system: Y (T*k x 1), X (T*k x k*p)
Y_vec = reshape(Y, (T - 1)* k, 1);
X_block = kron(eye(k), X);
Sigma_u_inv = inv(Sigma_u_hat);
W = kron(Sigma_u_inv, eye(T-1)); % Weighting matrix
beta_gls = (X_block' * W * X_block) \ (X_block' * W * Y_vec);


X=X_block;
Y=Y_vec;


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



    
    
    %% drawing Sigma and B|Sigma
    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*n,1) + reshape(PpsiTilde,n*m,1);
    Bdraw         = reshape(Bdraw,n*p+nex,n);
    % store reduced-form draws
    Bdraws{record,1}     = Bdraw;
    Sigmadraws{record,1} = Sigmadraw;