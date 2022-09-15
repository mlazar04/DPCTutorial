clear all;
close all;

%% DeePC data-driven predictive control software; authors P.C.N. Verheijen, M. Lazar 

%% Define Model (the user can define here any linear model)

[sysC, ymin, ymax, umin, umax, Mx, Nx, Mu, Nu] = PADC(); 

fs = 400e3;
Ts = 1/fs;

sys = c2d(sysC, Ts, 'zoh');
A = sys.A;
B = sys.B;
C = sys.C;
n = size(A,1);
ny = size(C,1);
nu = size(B,2);

%% computation of steady-state equilibrium values (only for plotting, i.e., they are not used in the cost function)
x_ss = Mx*5+Nx;
u_ss = Mu*5+Nu;

Wvar = 0.05; %noise variance for offline measured data
wvar = 0.05; %noise variance for online measured data

%% Controller parameters
N = 15;               % Prediction horizon
Q = 10;               % Output cost weight
R = 0.1*eye(2);       % Input cost weight
Tini = 12;            % Window of past input data used online for initializing the predicted trajectories




%% Generate data
% An offset / bias in tracking may still be observed if T is too low; this will also be observed for offset-free SPC with the same T
% The offset-free property is guaranteed only assuming the data lenth is large enough to ensure good estimation

T = 400; %the larger the data length, the better the robust performance; complexity of DeePC increases; 

%DeePC will take longer to simulate for the chosen T; T must be chosen large enough in
%relation to the noise variance to arrive at a reasonable bias, in general;
%reducing below 400 will increase the bias for 0.05 variance.


U = 0.8*idinput([T+N+Tini, nu], 'PRBS', [0, 1], [0, 1])';
l1=5e+4;               % Weight of g regularization (zero if there is no noise)
l2=1e+4;               %lambda 2 
l3=1e+4;               %lambda 3 

X = zeros(n, size(U,2));
Y = zeros(ny,size(U,2));
dU = zeros(nu, size(U,2));

%sensor noise
W = Wvar*randn(ny,size(U,2));

 for k=1:size(U,2)
     Y(:,k) = C*X(:,k) + W(:,k);
     if(k < size(U,2))
         X(:,k+1) = A*X(:,k) + B*U(:,k);
     end
 end

%% Generate delta U from U
dU(:,1) = U(:,1);
for i=2:size(dU,2)
    dU(:,i) = U(:,i)-U(:,i-1);
end

disp(['SNR output noise: ', num2str(snr(Y, W))]);

%% Define offset-free DeePC Controller

%Use dU instead of U to obtain the integral action system response
dUp = zeros(Tini*nu, T);
dUf = zeros(N*nu, T);
Yp = zeros(Tini*ny, T);
Yf = zeros(N*ny, T);

%Build Hankel matrices
for i = 1:Tini
    dUp((i-1)*nu+1:i*nu, :) = dU(:, i  :i+T-1);
    Yp((i-1)*ny+1:i*ny, :) =  Y(:, i+1:i+T);
end

for i = 1:N
    dUf((i-1)*nu+1:i*nu, :) = dU(:, i+Tini  :i+Tini+T-1);
    Yf((i-1)*ny+1:i*ny, :) =  Y(:, i+Tini+1:i+Tini+T);
end

%computing matrix weight for least squares regularization of g
PI = pinv([dUp; Yp; dUf])*[dUp; Yp; dUf];

%Build YALMIP offset-free DeePC problem
Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);

du = sdpvar(nu*N,1);
y = sdpvar(ny*N,1);
ref = sdpvar(ny*N,1);
g = sdpvar(T, 1);
sigma_du = sdpvar(nu*Tini,1);
sigma_y = sdpvar(ny*Tini,1);
du_ini = sdpvar(Tini*nu, 1);
y_ini = sdpvar(Tini*ny, 1);
u_km = sdpvar(nu,1); %previous input

% Build objective function
objective = (y-ref)'*Omega*(y-ref)+du'*Psi*du;    %tracking offset-free DeePC cost
objective = objective + l1*(g'*(eye(length(PI))-PI)'*(eye(length(PI))-PI)*g) + l2*(sigma_du'*sigma_du) + l3*(sigma_y'*sigma_y); %regularization terms

%above the term g'*(eye(length(PI))-PI)'*(eye(length(PI))-PI)*g can be
%replaced by g'*g to obtain the original DeePC regularization, which is
%however not consistent 

% Build constraints
constraints = [du_ini==dUp*g+sigma_du, y_ini==Yp*g+sigma_y, y==Yf*g, du==dUf*g]; %offset-free DeePC equality constraints
for k = 1:N 
    %standard input and output bounds, now with delta U in mind
    constraints = [constraints, ymin<=y(ny*(k-1)+1:ny*k)<=ymax, umin<=u_km+kron(ones(1,k), eye(nu))*du(1:nu*k)<=umax];
end

Parameters = {du_ini, y_ini, ref, u_km};
Outputs = {du, y};

%% Quadprog is a standard Matlab QP solver; Mosek is numerically more consistent but must be installed (academic license is free);
%% other solvers can be used with YALMIP
options = sdpsettings('solver', 'quadprog', 'verbose', 0, 'debug', 0);

controller = optimizer(constraints,objective,options, Parameters, Outputs);

%% Initialize Simulation
Tmax = 600*Ts;
t = 0:Ts:Tmax;
simLen = size(t,2);
x0 = zeros(5,1); %initial condition


%note that the reference sequence must be N steps longer to ensure the 
%last controller updates still have a future reference
r = zeros(ny, simLen+N);
r(:,5e-5/Ts:end) = 5;
d = zeros(ny, simLen);
d(:,6e-4/Ts:1e-3/Ts) = 1*-0.025; %active (multiply by 1) or inactive disturbance (multiply by 0)

y = zeros(ny, simLen);
u = zeros(nu, simLen);
du = zeros(nu, simLen);
x = zeros(n, simLen+1);
x(:,1) = x0;

w = wvar*randn(ny, simLen);

%% Simulation
nbytes = fprintf('time: 0 of %d', Tmax);
err = 0;

for k = 1:simLen
    
    while nbytes > 0
        fprintf('\b')
        nbytes = nbytes - 1;
    end
    nbytes = fprintf('processing %.3f of %d, QP status: %s', Ts*k, Tmax, yalmiperror(err));
    
    %Compute new (measured) output
    y(:,k) = C*x(:,k) + w(:,k);

    if(k >= Tini+1)
        %Update controller
        Rk = r(:, k+1:k+N);
        Rk = Rk(:);
        dU_ini = du(:, k-Tini:k-1);
        dU_ini = dU_ini(:); %flatten vector
        Y_ini = y(:,k-Tini+1:k);
        Y_ini = Y_ini(:);

        [Sol, err] = controller({dU_ini, Y_ini, Rk, u(:,max(1,k-1))});

        dUk = Sol{1};
        Yk = Sol{2};
        du(:,k) = dUk(1:nu);
    else
        %To prevent initial feasibility issues, control system in open loop
        % for the first Tini values.
        du(:,k) = 1*rand(nu, 1)-u(:,max(1, k-1)); 
    end
    
    %Update system
    u(:,k) = u(:,max(1,k-1)) + du(:,k);
       
    x(:,k+1) = A*x(:,k) + B*u(:,k) +[0;0;0;0;1]*d(:,k);
end

%% Display simulation results

figure();
ax1 = subplot(311);
stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference');
hold on;
stairs(t, y(1,:), 'r','LineWidth',1, 'DisplayName', 'iDeePC');

ylabel('Output current [A]');
legend;

axis([0 t(end) 0 5.5]);

ax2 = subplot(312);
stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference');
hold on;
stairs(t, y(1,:), 'r', 'LineWidth',1, 'DisplayName', 'iDeePC');

ylabel('Output current [A]');
axis([0 t(end) 4.5 5.5]);

ax3 = subplot(313);
hold on;
stairs(t, u(1,:), 'r','LineWidth',1);
stairs(t, u_ss(1)*ones(1,length(t)), '--k');
stairs(t, u(2,:), 'b','LineWidth',1);
stairs(t, u_ss(2)*ones(1,length(t)), '--k');

ylabel('Duty-cycles ');
xlabel('Time [s]');

