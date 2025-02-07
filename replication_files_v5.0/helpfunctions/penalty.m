function y= penalty(q)

global ssigma;
global objective;



x2 = -(1/ssigma(2,1))*objective*q;

if x2>0
    y2=100*x2;
else
    y2=x2;
end

y=y2;


