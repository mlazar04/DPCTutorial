%provides a time scale plot of the Hanke-Raus heuristic for tuning
%regularized iDeePC (alternative to the function DeePCRegularizationPlot)
%requires the N, Q, R, ry, ru, and the data matrices dUp, dUf, Yp, Yf
%can also be used for tuning regularized DeePC by using the matrices Up, Uf
%instead of dUp, dUf, respectively
%authors: M. Lazar, P.C.N. Verheijen

fh=[];
lamh=[];
N=15;
Q = 10;               % Output cost weight
R = 0.1*eye(2);  
Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);

%the data matrices come after runing iDeePC code
At=[chol(Omega)*Yf; chol(Psi)*dUf];
PI = pinv([dUp; Yp; dUf])*[dUp; Yp; dUf];

ry=5*ones(N,1);
ru=0*ones(2*N,1);

bt=[chol(Omega)*ry; chol(Psi)*ru];

L=eye(length(PI))-PI;
Linv=pinv(L);
At=At*Linv;

%Hanke and Raus heuristic     
for lam=10:10:100000
    aux=inv(At*At'+lam*eye(length(At*At')));
    flam=sqrt(bt'*lam^2*aux^3*bt);
    lamh=[lamh lam];
    fh=[fh flam];
end
plot(lamh,fh,'-ob');
hold on;

%Tikhonov heuristic
lamh=[];
fh=[];
for lam=10:10:100000
    aux=inv(At*At'+lam*eye(length(At*At')));
    flam=sqrt(bt'*lam^2*aux^4*At*At'*bt);
    lamh=[lamh lam];
    fh=[fh flam];
    end    


plot(lamh,fh,'-or');

legend('Hanke-Raus', 'Tikhonov');