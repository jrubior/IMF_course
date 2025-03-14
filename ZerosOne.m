close all
clear
clc
T = 1000;
n_draws=100;

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
X=[ones(T-1,1) Y(1:end-1,:)]; % this is the notation used in our papers
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
Adraws=cell([n_draws,2]);
Ldraws=cell([n_draws,p+1]);

accdraws=cell([n_draws,3]);
accAdraws=cell([n_draws,2]);
accLdraws=cell([n_draws,p+1]);


%% Signs

S=zeros(1,n);
S(1,1)=1; % first variable
ee=eye(n);
shock=ee(:,1); % first shock

Z=zeros(1,n);
Z(1,2)=1; % second variable

for i=1:n_draws

    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*n,1) + reshape(PpsiTilde,n*m,1);
    Bdraw         = reshape(Bdraw,n*p+nex,n);
    draws{i,1}    = Bdraw;
    draws{i,2}    = Sigmadraw;
    %draws{i,3}    = DrawQ(n);

    Q = zeros(n, n);

    j=1;

    nz=size(Z,1);

    x_j = randn(n + 1 - j-nz, 1);

    % Normalize to get w_j
    w_j = x_j / norm(x_j);

    M_j = Z*hh(draws{i,2})';
    K_j= null(M_j);

    % Compute q_j
    q_j = K_j * w_j;

    % Store q_j in Q
    Q(:, j) = q_j;


    j=2;
    nz=size(Z,1);

    x_j = randn(n + 1 - j, 1);

    % Normalize to get w_j
    w_j = x_j / norm(x_j);

    M_j = Q(:, 1:(j-1))';
    K_j= null(M_j);

    % Compute q_j
    q_j = K_j * w_j;

    % Store q_j in Q
    Q(:, j) = q_j;


    draws{i,3}= Q;

    Adraws{i,1}=hh(draws{i,2})\draws{i,3};
    Adraws{i,2}=draws{i,1}*Adraws{i,1};

    Ldraws{i,1}=hh(draws{i,2})'*draws{i,3};
    Ldraws{i,2}=draws{i,1}(nex+1:nex+2,:)'*Ldraws{i,1};


    if S*Ldraws{i,1}*shock > 0

        accdraws{i,1}    = draws{i,1};
        accdraws{i,2}    = draws{i,2};
        accdraws{i,3}    = draws{i,3};

        accAdraws{i,1}=hh(accdraws{i,2})\accdraws{i,3};
        accAdraws{i,2}=accdraws{i,1}*accAdraws{i,1};

        accLdraws{i,1}=hh(accdraws{i,2})'*accdraws{i,3};
        accLdraws{i,2}=accdraws{i,1}(nex+1:nex+2,:)'*accLdraws{i,1};

    end

end

% Remove empty cells from accdraws, accAdraws, and accLdraws

accdraws = accdraws(~cellfun(@isempty, accdraws(:,1)), :);
accAdraws = accAdraws(~cellfun(@isempty, accAdraws(:,1)), :);
accLdraws = accLdraws(~cellfun(@isempty, accLdraws(:,1)), :);

disp(['How many do we keep? ', num2str(size(accLdraws,1))]);

