close all
clear
clc

load datadata.mat

traindimvalue=0.8;

nxx=size(xx,1);
traindim=round(nxx*traindimvalue);
checkdim=nxx-traindim;

% Create a pattern recognition network
hiddenLayerSize = 20;
net = feedforwardnet(hiddenLayerSize);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 0.8; % 70% for training
net.divideParam.valRatio = 0.1;  % 15% for validation
net.divideParam.testRatio = 0.1; % 15% for testing

xx=xx';
yy=yy';
yy1=yy1';
yy2=yy2';

xxtrain=xx(:,1:traindim);
yytrain=yy(:,1:traindim);
yy1train=yy1(:,1:traindim);
yy2train=log(yy2(:,1:traindim));

% Train the Network

[net,tr] = train(net,xxtrain,yy2train,'useParallel','yes');

xxtest=xx(:,traindim+1:end);
yytest=yy(:,traindim+1:end);
yy1test=yy1(:,traindim+1:end);
yy2test=log(yy2(:,traindim+1:end));

yy2net=zeros(1,size(xxtest,2));

for i=1:size(xxtest,2)

    yy2net(1,i)=net(xxtest(:,i));

end

toc

close all
plot(1:size(xxtest,2),yy2test,'color','r');
hold on
plot(1:size(xxtest,2),yy2net,'color','b');


corr(yy2test',yy2net')



