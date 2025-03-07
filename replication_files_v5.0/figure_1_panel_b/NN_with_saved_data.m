close all
clear
clc

load datadata.mat

traindimvalue=0.9;

nxx=size(xx,1);
traindim=round(nxx*traindimvalue);
checkdim=nxx-traindim;

% Create a pattern recognition network
hiddenLayerSize = 10;
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

uwtest=exp(yy2test)/sum(exp(yy2test));
uwnet=exp(yy2net)/sum(exp(yy2net));

close all
plot(1:size(xxtest,2),uwtest,'color','r');
hold on
plot(1:size(xxtest,2),uwnet,'color','b');

corr(uwtest',uwnet')
cov(uwtest,uwnet)

save net.mat net xxtest yytest yy1test yy2test xxtrain yytrain yy1train yy2train 



