function DeePCRegularizationPlot(Up, Yp, Uf, Yf, Omega, Psi, Ry, Ru, range)
    %DeePCRegularizationPlot(Up, Yp, Uf, Yf, Omega, Psi, Ry, Ru, range)
    %Authors: P.C.N. Verheijen, M. Lazar
    % plots the regularization score for various values
    % of lambda. Low scores for lambda push the estimate towards the
    % unbiased estimate. In general, any lambda << min(Psi, Omega) should
    % be disregarded.
    %
    % Usage in DeePC problem:
    % min_{U,Y,g} (Y-Ry)'*Omega*(Y-Ry) + (U-Ru)'*Psi*(U-Ru)
    %             + lambda*norm((I-PI)*g)
    % s.t. [Up; Yp; Uf; Yf]*g = [u_ini; y_ini; U; Y]
    %      {inequality constraints}
    %
    % where PI = pinv([Up; Yp; Uf])*[Up; Yp; Uf];
    % the last argument, range = [min, max, nPts], can be used to specify the
    % range of the plot. This function generates 2 plots, 1 with
    % logarithmic spacing, one with linear spacing. Each having nPts points
    % between min and max. If range is omitted, 
    % range = [1e-1*min(Psi, Omega), 1e4*max(Psi, Omega), 1000]

    if(nargin == 8)
        range = [1e-1*min(min(Psi(Psi>0)), min(Omega(Omega>0))), 1e4*max(max(diag(Psi)), max(diag(Omega))), 1000];
    end
    
    %the data matrices come after running iDeePC code
    At=[chol(Omega)*Yf; chol(Psi)*Uf];
    PI = pinv([Up; Yp; Uf])*[Up; Yp; Uf];

    bt=[chol(Omega)*Ry; chol(Psi)*Ru];

    L=eye(length(PI))-PI;
    
    % Use a higher tolerance here, after all, the eigenvalues of L are
    % either 1 or 0.
    Linv=pinv(L, 1e-8);
    At=At*Linv;

    %% First evaluate the heuristics using a logarithmic spacing
    fh=[];
    lamLog = logspace(log10(range(1)), log10(range(2)), range(3));
    %Hanke and Rause heuristic     
    for lam=lamLog
        aux=inv(At*At'+lam*eye(length(At*At')));
        flam=sqrt(bt'*lam^2*aux^3*bt);
        fh=[fh flam];
    end
    
    %Tikhonov heuristic
    fhT=[];
    for lam=lamLog
        aux=inv(At*At'+lam*eye(length(At*At')));
        flam=sqrt(bt'*lam^2*aux^4*At*At'*bt);
        fhT=[fhT flam];
    end
    figure();
    semilogx(lamLog,fh,'b', 'LineWidth', 2, 'DisplayName', 'Hanke and Rause');
    hold on;
    semilogx(lamLog,fhT,'r', 'LineWidth', 2, 'DisplayName', 'Tikhonov');
    title('Logarithmic regularization score plot');
    ylabel('Regularization score');
    xlabel('lambda');
    grid on;
    legend;
    
    %% Now the same for a linear spacing
    fh=[];
    lamLin = linspace(range(1), range(2), range(3));
    %Hanke and Rause heuristic     
    for lam=lamLin
        aux=inv(At*At'+lam*eye(length(At*At')));
        flam=sqrt(bt'*lam^2*aux^3*bt);
        fh=[fh flam];
    end

    %Tikhonov heuristic
    fhT=[];
    for lam=lamLin
        aux=inv(At*At'+lam*eye(length(At*At')));
        flam=sqrt(bt'*lam^2*aux^4*At*At'*bt);
        fhT=[fhT flam];
    end
    figure();
    plot(lamLin,fh,'b', 'LineWidth', 2, 'DisplayName', 'Hanke and Rause');
    hold on;
    plot(lamLin,fhT,'r', 'LineWidth', 2, 'DisplayName', 'Tikhonov');
    ylabel('Regularization score');
    xlabel('lambda');
    title('Linear regularization score plot');
    grid on;
    legend;
end