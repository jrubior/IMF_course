%==========================================================================
%% housekeeping
%==========================================================================
clear variables;
close all;
userpath('clear');
clc;
tic;

rng('default'); % reinitialize the random number generator to its startup configuration
rng(0);         % set seed

message = 'Please wait while Panel (b) of Figure (1) is being replicated ... ';
disp(message);

currdir='/Users/jrubior/Documents/GitHub/IMF_course/replication_files_v5.0/figure_1_panel_b';
cd ..
get_help_dir_currdir='/Users/jrubior/Documents/GitHub/IMF_course/replication_files_v5.0';
addpath([get_help_dir_currdir,'/helpfunctions']); % set path to helper functions
cd(currdir)


%==========================================================================
%% load the data and priors
%==========================================================================
% variables are in the following order:
% 1: Adjusted TFP
% 2: Stock Prices
% 3: Consumption
% 4: Real Interest Rate
% 5: Hours Worked
data = importdata([get_help_dir_currdir,'/data/data.csv']); % import Beaudry, Nam and Wang (2011) data
num  = data.data(:,3:7)*100;        % get variables to be used in SVAR

%==========================================================================
%% model setup
%==========================================================================
nlag      = 4;               % number of lags
nvar      = 5;               % number of endogenous variables
nex       = 1;               % set equal to 1 if a constant is included; 0 otherwise
m         = nvar*nlag + nex; % number of exogenous variables
nd        = 1e5;             % number of orthogonal-reduced-form (B,Sigma,Q) draws
iter_show = 1e3;            % display iteration every iter_show draws
horizon   = 40;              % maximum horizon for IRFs
index     = 40;              % define  horizons for the FEVD
NS        = 1;               % number of objects in F(THETA) to which we impose sign and zero restrictios: F(THETA)=[L_{0}]
e         = eye(nvar);       % create identity matrix
maxdraws  = 1e2;             % max number of importance sampling draws
conjugate = 'structural';    % structural or irfs or empty

%==========================================================================
%% identification: declare Ss and Zs matrices
%==========================================================================
% restrictions on IRFs
% horizons to restrict
horizons = 0;

% sign restrictions
S = cell(nvar,1);
for ii=1:nvar
    S{ii}=zeros(0,nvar*NS);
end
ns1  = 1; % ns_{1}  = 1  positive sign restriction on the contemporanouos response of Stock Prices (variable 2) to a unit optimism shock
S{1} = zeros(ns1,nvar*NS);
S{1}(1,2) = 1;



% zero restrictions
Z=cell(nvar,1);
for i=1:nvar
    Z{i}=zeros(0,numel(  horizons)*nvar);
end

nz1=1; % nz_{1}  = 1 zero restriction on the contemporanouos response of TFP (variable 1) to a unit optimism shock
Z{1}=zeros(nz1,numel(horizons)*nvar);
Z{1}(1,1)=1;


%==========================================================================
%% Setup info
%==========================================================================
info=SetupInfo(nvar,m,Z,@(x)chol(x));

% ZIRF()
info.nlag     = nlag;
info.horizons = horizons;
info.ZF       = @(x,y)ZIRF(x,y);

% functions useful to compute the importance sampler weights
fs      = @(x)ff_h(x,info);
r       = @(x)ZeroRestrictions(x,info);

if strcmp(conjugate,'irfs')==1
    fo              = @(x)f_h(x,info);
    fo_str2irfs     = @(x)StructuralToIRF(x,info);
    fo_str2irfs_inv = @(x)IRFToStructural(x,info);
    r_irfs          = @(x)IRFRestrictions(x,info); 
end


% function useful to check the sign restrictions
fh_S_restrictions  = @(y)StructuralRestrictions(y,S);

%==========================================================================
%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [yt(t-1), ... , yt(t-nlag), constant];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
yt = num(nlag+1:end,:);
T  = size(yt,1);
xt = zeros(T,nvar*nlag+nex);
for i=1:nlag
    xt(:,nvar*(i-1)+1:nvar*i) = num((nlag-(i-1)):end-i,:) ;
end
if nex>=1
    xt(:,nvar*nlag+nex)=ones(T,1);
end
% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors


%% prior for reduced-form parameters
nnuBar              = 0;
OomegaBarInverse    = zeros(m);
PpsiBar             = zeros(m,nvar);
PphiBar             = zeros(nvar);

%% posterior for reduced-form parameters
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
PpsiTilde           = OomegaTilde*(X'*Y + OomegaBarInverse*PpsiBar);
PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - PpsiTilde'*OomegaTildeInverse*PpsiTilde;
PphiTilde           = (PphiTilde'+PphiTilde)*0.5;


%% useful definitions
% definitios used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
Bdraws         = cell([nd,1]); % reduced-form lag parameters
Sigmadraws     = cell([nd,1]); % reduced-form covariance matrices
Qdraws         = cell([nd,1]); % orthogonal matrices
storevefh      = zeros(nd,1);  % volume element f_{h}
storevegfhZ    = zeros(nd,1);  % volume element g o f_{h}|Z
uw             = zeros(nd,1);  % unnormalized importance sampler weights
structuraldraws  = zeros(nd,nvar*nvar+nvar*m); % structural lag parameters

if strcmp(conjugate,'irfs')==1
    storevephi      = zeros(nd,1);  % volume element f_{h}
    storevegphiZ    = zeros(nd,1);  % volume element g o f_{h}|Z
end

% definitions related to IRFs; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
J      = [e;repmat(zeros(nvar),nlag-1,1)];
A      = cell(nlag,1);
extraF = repmat(zeros(nvar),1,nlag-1);
F      = zeros(nlag*nvar,nlag*nvar);
for l=1:nlag-1
    F((l-1)*nvar+1:l*nvar,nvar+1:nlag*nvar)=[repmat(zeros(nvar),1,l-1) e repmat(zeros(nvar),1,nlag-(l+1))];
end

% definition to facilitate the draws from B|Sigma
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

%% initialize counters to track the state of the computations
counter = 1;
record  = 1;
count   = 0;

tic

while record<=nd


    
    
    %% step 1 in Algorithm 2
    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(PpsiTilde,nvar*m,1);
    Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
    % store reduced-form draws
    Bdraws{record,1}     = Bdraw;
    Sigmadraws{record,1} = Sigmadraw;
    
   
    %% steps 2:4 of Algorithm 2
    w           = DrawW(info);   
    x           = [vec(Bdraw); vec(Sigmadraw); w];
    structpara  = ff_h_inv(x,info);
    structuraldraws(record,:) = structpara;
    
    % store the matrix Q associated with step 3
    Qdraw            = SpheresToQ(w,info,Bdraw,Sigmadraw);
    Qdraws{record,1} = reshape(Qdraw,nvar,nvar);


    %% check if sign restrictions hold
    signs      = fh_S_restrictions(structpara);
    
    
    if (sum(signs>0))==size(signs,1)
        
        count=count+1;
  
        %% compute importance sampling weights
        
        switch conjugate
            
            case 'structural'
                
                
                storevefh(record,1)   = (nvar*(nvar+1)/2)*log(2)-(2*nvar+m+1)*LogAbsDet(reshape(structpara(1:nvar*nvar),nvar,nvar));
                storevegfhZ(record,1) = LogVolumeElement(fs,structpara,r); 
                uw(record,1)          = exp(storevefh(record,1) - storevegfhZ(record,1));
                
            case 'irfs'
                
                irfpara                = fo_str2irfs(structpara);
                storevephi(record,1)   = LogVolumeElement(fo,structpara)   + LogVolumeElement(fo_str2irfs_inv,irfpara);
                storevegphiZ(record,1) = LogVolumeElement(fs,structpara,r) + LogVolumeElement(fo_str2irfs_inv,irfpara,r_irfs); 
                uw(record,1)           = exp(storevephi(record,1) - storevegphiZ(record,1));
                
            otherwise
                
                uw(record,1) = 1;
                
        end
        
    else
        
        uw(record,1) = 0;
        
    end
    
    if counter==iter_show

        toc
        tic 
        
        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(nd-(record))])
        counter =0;
        
    end
    counter = counter + 1;
    
    record=record+1;
    
end
toc


% compute the normalized weights and estimate the effective sample size of the importance sampler
imp_w  = uw/sum(uw);
ne     = floor(1/sum(imp_w.^2));


%% useful definitions to store relevant objects
A0tilde       = zeros(nvar,nvar,ne);               % define array to store A0
Aplustilde    = zeros(m,nvar,ne);                  % define array to store Aplus
Ltilde        = zeros(horizon+1,nvar,nvar,ne);     % define array to store IRF
FEVD          = zeros(nvar,ne);                 % define array to store FEVD
% initialize counter to track the state of the importance sampler
count_IRF     = 0;
for s=1:min(ne,maxdraws)
    
    %% draw: B,Sigma,Q
    is_draw     = randsample(1:size(imp_w,1),1,true,imp_w);
    Bdraw       = Bdraws{is_draw,1};
    Sigmadraw   = Sigmadraws{is_draw,1};
    Qdraw       = Qdraws{is_draw,1};
    
    x          = [reshape(Bdraw,m*nvar,1); reshape(Sigmadraw,nvar*nvar,1); Qdraw(:)];
    structpara = f_h_inv(x,info);
    

    LIRF = IRF_horizons(structpara, nvar, nlag, m, 0:horizon);
    
    for h=0:horizon
        Ltilde(h+1,:,:,s) =  LIRF(1+h*nvar:(h+1)*nvar,:);
    end
    
    % FEVD
    % compute matrix F: useful for FEVD
    hSigmadraw = hh(Sigmadraw);
    A0         = hSigmadraw\e;
    Aplus      = Bdraw*A0;
    % Page 8 ARRW
    for l=1:nlag-1
        A{l} = Aplus((l-1)*nvar+1:l*nvar,1:end);
        F((l-1)*nvar+1:l*nvar,1:nvar)=A{l}/A0;
    end
    A{nlag} = Aplus((nlag-1)*nvar+1:nlag*nvar,1:end);
    F((nlag-1)*nvar+1:nlag*nvar,:)=[A{nlag}/A0 extraF];
    
    FEVD(:,s)         = variancedecomposition(F',J,Sigmadraw, hh(Sigmadraw)'*Qdraw(:,1),nvar,index);
    
    % store weighted independent draws
    A0tilde(:,:,s)    = reshape(structpara(1:nvar*nvar),nvar,nvar);
    Aplustilde(:,:,s) = reshape(structpara(nvar*nvar+1:end),m,nvar);
    
end
A0tilde    = A0tilde(:,:,1:s);
Aplustilde = Aplustilde(:,:,1:s);
Ltilde     = Ltilde(:,:,:,1:s);
FEVD       = FEVD(:,1:s);

addpath('results/plothelpfunctions'); 
store_results_and_plot_IRFs;

message = 'Done.';
disp(message);

%% Neural Network Toolbox

nonzero_idx = storevefh ~= 0;

% Keep only the rows where y is nonzero
y = storevefh(nonzero_idx, :);
x = structuraldraws(nonzero_idx, :);

nx=size(x,2);

% Define Neural Network Architecture
layers = [
    featureInputLayer(nx)                    % Input layer (1 feature)
    fullyConnectedLayer(10)                  % Hidden layer with 10 neurons                 % Another hidden layer
             % Hidden layer with 10 neurons                 % Another hidden layer
    reluLayer                             % Activation function
    fullyConnectedLayer(1)                    % Output layer (1 neuron)
    regressionLayer                           % Regression output
];

initialLR = 0.1;  % Initial learning rate
lambda = 0.01;    % Decay rate
numEpochs = 50;   % Total training epochs
miniBatchSize = 2^11;


% Specify Training Options
 options = trainingOptions('adam', ...
     'MaxEpochs', numEpochs, ...
     'MiniBatchSize', miniBatchSize, ...
     'InitialLearnRate', initialLR, ...
     'Verbose', true);

% Train the Network
net = trainNetwork(x, y, layers, options);

save('net.mat', 'net');





