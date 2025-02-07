function output = variancedecomposition(A,J,Ssigmau,P,n,t)
    MSE = zeros(n,n);
    for i=0:(t-1)
        Pphi = (J'*(A^i)*J);
        MSE  = Pphi*Ssigmau*Pphi' + MSE;
    end
    MSEdiag  = diag(MSE);
    W   = zeros(n,1);
    evp = zeros(n,n);
    for i=0:(t-1)
        Pphi   = (J'*(A^i)*J);
        Ttheta = Pphi*P;
        evp    =  Ttheta*Ttheta' + evp;
    end
    W(:,1) = diag(evp)./MSEdiag;
    output = W(:,1)';
end
