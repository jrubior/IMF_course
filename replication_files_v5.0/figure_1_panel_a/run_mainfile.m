%% housekeeping
clear variables;
close all;
userpath('clear');
clc;
tic;
global ssigma objective;


rng('default'); % reinitialize the random number generator to its startup configuration
rng(0);         % set seed


message = 'Please wait while Panel (a) of Figure 1 is being replicated ... ';
disp(message);


currdir=pwd;
cd ..
get_help_dir_currdir=pwd;
addpath([get_help_dir_currdir,'/helpfunctions']); % set path to helper functions
cd(currdir)



%% load the data and priors
% variables are in the following order:
% 1: Adjusted TFP
% 2: Stock Prices
% 3: Consumption
% 4: Real Interest Rate
% 5: Hours Worked
data = importdata([get_help_dir_currdir,'/data/data.csv']); % import Beaudry, Nam and Wang (2011) data
num  = data.data(:,3:7)*100;        % get variables to be used in SVAR


%% model setup
p   = 4;            % number of lags
n   = 5;            % number of endogenous variables
nex = 1;            % set equal to 1 if a constant is included; 0 otherwise
m   = n*p + nex;    % number of independent variables
nd  = 1e4;          % number of orthogonal-reduced-form (B,Sigma,Q) draws
iter_show = nd/5e0; % display iteration every iter_show draws
horizon = 40;       % maximum horizon for IRFs
index   = 40;       % define  horizons for the FEVD
NS = 1;             % number of objects in F(THETA) to which we impose sign and zero restrictios: F(THETA)=[L_{0}]
e  = eye(n);        % create identity matrix
maxdraws = 1e4;     % max number of importance sampling draws
agnostic = 'irfs';  % select: 'irfs' or structural;


%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [yt(t-1), ... , yt(t-p), constant];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
yt = num(p+1:end,:);
T  = size(yt,1);
xt = zeros(T,n*p+nex);
for i=1:p
    xt(:,n*(i-1)+1:n*i) = num((p-(i-1)):end-i,:) ;
end
if nex>=1
    xt(:,n*p+nex)=ones(T,1);
end
% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by n matrix of observations
X = xt; % T by (n*p+1) matrix of regressors
B = (X'*X)\(X'*Y); % ols estimates
U = Y-X*B;         % ols residuals
Sigmau = U'*U/T;   % ols covariance matrix of residuals
ssigma = sqrt(diag(Sigmau)); % scale for penalty function 


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


%% set help functions
hh       = @(x)chol(x);


%% useful definitions
% definitios used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
Bdraws       = cell([nd,1]); % reduced-form lag parameters
Sigmadraws   = cell([nd,1]); % reduced-form covariance matrices
Qdraws       = cell([nd,1]); % orthogonal matrices
Ltilde       = zeros(horizon+1,n,nd);     % define array to store IRF
FEVD         = zeros(n,nd); % define array to store FEVD

% definitions related to IRFs; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
J      = [e;repmat(zeros(n),p-1,1)];
A      = cell(p,1);
extraF = repmat(zeros(n),1,p-1);
F      = zeros(p*n,p*n);
for l=1:p-1
    F((l-1)*n+1:l*n,n+1:p*n) = [repmat(zeros(n),1,l-1) e repmat(zeros(n),1,p-(l+1))];
end

% definition to facilitate the draws from B|Sigma
cholOomegaTilde=hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

%% initialize counters to track the state of the computations
counter = 1;
record  = 1;


while record<=nd
    
    
    %% drawing Sigma and B|Sigma
    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*n,1) + reshape(PpsiTilde,n*m,1);
    Bdraw         = reshape(Bdraw,n*p+nex,n);
    % store reduced-form draws
    Bdraws{record,1}     = Bdraw;
    Sigmadraws{record,1} = Sigmadraw;
    
    
    %% Mountford and Uhlig (2009)'s PFA
    objective = e(2,:)*hh(Sigmadraw)';
    Aeq       = e(1,:)*hh(Sigmadraw)';
    beq       = 0;
    q1ga      = rand(n,1);
    [q,valq]  = fmincon(@penalty,q1ga,[],[],Aeq,beq,[],[],@mycon,optimset('MaxFunEvals',40000,'MaxIter',20000,'Display','off','Algorithm','active-set'));
    % compute matrix F: useful for computing IRFs and FEVD
    hSigmadraw = hh(Sigmadraw);
    A0         = hSigmadraw\e;
    Aplus      = Bdraw*A0;
    % definitions related to IRFs; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
    for l=1:p-1
        A{l} = Aplus((l-1)*n+1:l*n,1:end);
        F((l-1)*n+1:l*n,1:n)=A{l}/A0;
    end
    A{p} = Aplus((p-1)*n+1:p*n,1:end);
    F((p-1)*n+1:p*n,:)=[A{p}/A0 extraF];
    for h=1:horizon+1
        Ltilde(h,:,record) = (J'*((F')^(h-1))*J)*hSigmadraw'*q;
    end
    FEVD(:,record) = variancedecomposition(F',J,Sigmadraw,hh(Sigmadraw)'*q,n,index);

    if counter==iter_show
        
        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(nd-(record))])
        counter =0;
        
    end
    counter = counter + 1;
    record = record +1;
end

addpath('results/plothelpfunctions'); 
store_results_and_plot_IRFs;

message = 'Done.';
disp(message);
