close all
clear
clc

load datadata.mat

traindimvalue=0.8;

nxx=size(xx,1);
traindim=round(nxx*traindimvalue);
checkdim=nxx-traindim;

% Create a pattern recognition network
hiddenLayerSize = [10 10];
net = feedforwardnet(hiddenLayerSize);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 0.7; % 70% for training
net.divideParam.valRatio = 0.15;  % 15% for validation
net.divideParam.testRatio = 0.15; % 15% for testing

xx=xx';
yy=yy';
yy1=yy1';
yy2=yy2';

xxtrain=xx(:,1:traindim);
yytrain=yy(:,1:traindim);
yy1train=yy1(:,1:traindim);
yy2train=yy2(:,1:traindim);

% Train the Network

[net,tr] = train(net,xxtrain,yytrain,'useParallel','yes');

xxtest=xx(:,traindim+1:end);
yytest=yy(:,traindim+1:end);
yy1test=yy1(:,traindim+1:end);
yy2test=yy2(:,traindim+1:end);

yynet=zeros(1,size(xxtest,2));

for i=1:size(xxtest,2)

    yynet(1,i)=net(xxtest(:,i));

end

errornet=abs(yytest-yynet);
toc

close all
plot(1:size(xxtest,2),yytest,'color','r');
hold on
plot(1:size(xxtest,2),yynet,'color','b');


corr(ytest',yynet')



