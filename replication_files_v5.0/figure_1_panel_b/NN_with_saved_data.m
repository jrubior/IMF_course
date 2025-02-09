close all
clear
clc

load net.mat

nx=size(x,2);

% % Define Neural Network Architecture
% layers = [
%     featureInputLayer(nx)                    % Input layer (1 feature)
%     fullyConnectedLayer(20)                  % Hidden layer with 10 neurons                 % Another hidden layer
%              % Hidden layer with 10 neurons                 % Another hidden layer
%     sigmoidLayer                             % Activation function
%     fullyConnectedLayer(10)                    % Output layer (1 neuron)
%     sigmoidLayer
%     fullyConnectedLayer(1)
%     regressionLayer% Regression output
% ];
%
% initialLR = 0.1;  % Initial learning rate
% lambda = 0.01;    % Decay rate
% numEpochs = 50;   % Total training epochs
% miniBatchSize = 2^13;
%
%
% % Specify Training Options
%  options = trainingOptions('adam', ...
%      'MaxEpochs', numEpochs, ...
%      'MiniBatchSize', miniBatchSize, ...
%      'InitialLearnRate', initialLR, ...
%      'LearnRateSchedule', 'piecewise', ...
%      'Shuffle', 'every-epoch', ...
%      'Verbose', true);
%
% % Train the Network
% net = trainNetwork(x, y, layers, options);
%
% save('net.mat', 'net', 'x', 'y');
%
%
% % Retrain the Network
% [net, info] = trainNetwork(x, y, net.Layers, options);
% loss_save=info.TrainingRMSE(end)

% % Define Neural Network Architecture
% layers = [
%     featureInputLayer(nx)   % Input layer with nx features (adjust automatically)
%     fullyConnectedLayer(20)          % First hidden layer with 20 neurons
%     sigmoidLayer                     % Log-Sigmoid activation (logsig equivalent)
%     fullyConnectedLayer(10)          % Second hidden layer with 10 neurons
%     sigmoidLayer                     % Log-Sigmoid activation (logsig equivalent)
%     fullyConnectedLayer(1)           % Output layer (1 neuron, behaves as purelin)
%     regressionLayer                   % Regression output (purelin behavior)
% ];
%
% % Training Options (Equivalent to trainRatio, valRatio, testRatio)
% options = trainingOptions('adam', ...  % Equivalent to feedforwardnet default
%     'MaxEpochs', 1000, ...                % Default epoch limit in feedforwardnet
%     'ValidationPatience', 6, ...           % Equivalent to early stopping
%     'InitialLearnRate', 0.01, ...          % Typical learning rate for LM
%     'Shuffle', 'every-epoch', ...
%     'Verbose', true, ...
%     'Plots', 'training-progress');
%
% % Train the Network
% net = trainNetwork(x, y, layers, options);

% Create a pattern recognition network
hiddenLayerSize = [20 10];

net = feedforwardnet(hiddenLayerSize);
%net.layers{2}.transferFcn = 'purelin';


% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 0.7; % 70% for training
net.divideParam.valRatio = 0.15;  % 15% for validation
net.divideParam.testRatio = 0.15; % 15% for testing

x=x(:,1:25);

xtest=x(40001:end,:);
ytest=y(40001:end,:);
x=x(1:40000,:);
y=y(1:40000,:);

x=x';
y=y';

% Train the Network

[net,tr] = train(net,x,y);

xtest=xtest';
ytest=ytest';

ynet=zeros(1,size(xtest,2));

for i=1:size(xtest,2)

    ynet(1,i)=net(xtest(:,i));

end

errornet=abs(ytest-ynet);
toc

close all
plot(1:size(xtest,2),ytest,'color','r');
hold on
plot(1:size(xtest,2),ynet,'color','b');


corr(ytest',ynet')



