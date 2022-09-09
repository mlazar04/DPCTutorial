clear all;
close all;

%% Offset-free SPC data-driven predictive control software; authors P.C.N. Verheijen, M. Lazar 

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

Wvar = 0.1; %noise variance for offline measured data
wvar = 0.1; %noise variance for online measured data

%% Controller parameters
N = 15;               % Prediction horizon
Q = 10;               % Output cost weight
R = 1*eye(2);       % Input cost weight

Tini = 12;    % Window of past input data used online for initializing the predicted trajectories




%% Generate data
T = 4000;  %size of the data sequence used for estimating the SPC prediction matrices; complexity of the SPC algorithm is not affected;
U = 0.8*idinput([T+N+Tini, nu], 'PRBS', [0, 1], [0, 1])';
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

%% Define offset-free SPC Controller

%Use dU instead of U to obtain the integral action matrices
dUp = zeros(Tini*nu, T);
dUf = zeros(N*nu, T);
Yp = zeros(Tini*ny, T);
Yf = zeros(N*ny, T);


for i = 1:Tini
    dUp((i-1)*nu+1:i*nu, :) = dU(:, i  :i+T-1);
    Yp((i-1)*ny+1:i*ny, :) = Y(:, i+1:i+T);
end

for i = 1:N
    dUf((i-1)*nu+1:i*nu, :) = dU(:, i+Tini  :i+Tini+T-1);
    Yf((i-1)*ny+1:i*ny, :) = Y(:, i+Tini+1:i+Tini+T);
end


Theta = Yf*pinv([dUp;Yp;dUf]);
P1    = Theta(:, 1:Tini*nu);
P2    = Theta(:, Tini*nu+1:Tini*(nu+ny));
Gamma = Theta(:, Tini*(nu+ny)+1:end);

%Build YALMIP offset-free SPC problem
Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);

du = sdpvar(nu*N,1);
y = sdpvar(ny*N,1);
ref = sdpvar(ny*N,1);
du_ini = sdpvar(Tini*nu, 1);
y_ini = sdpvar(Tini*ny, 1);
u_km = sdpvar(nu,1); %previous input
objective = (y-ref)'*Omega*(y-ref)+du'*Psi*du;
constraints = [y == P1*du_ini+P2*y_ini+Gamma*du];
for k = 1:N
    constraints = [constraints, ymin<=y(ny*(k-1)+1:ny*k,:)<=ymax;  umin<=u_km+kron(ones(1,k), eye(nu))*du(1:nu*k)<=umax];
end

Parameters = {du_ini, y_ini, ref, u_km};
Outputs = {du, y};

%% Quadprog is a standard Matlab QP solver; Mosek is numerically more consistent; other solvers can be used with YALMIP
options = sdpsettings('solver', 'quadprog', 'verbose', 0, 'debug', 0);

controller = optimizer(constraints,objective,options, Parameters, Outputs);

%% Initialize Simulation (the simulation duration and plotting part should be adapted to each example / system)
Tmax = 600*Ts;
t = 0:Ts:Tmax;
simLen = size(t,2);
x0 = zeros(5,1); %initial condition

%reference and disturbance sequences. Note that reference must be N samples
%longer to accommodate for the last predictions
r = zeros(ny, simLen+N);
r(:,5e-5/Ts:end) = 5;
d = zeros(ny, simLen);
d(:,6e-4/Ts:1e-3/Ts) = -0.025;

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
    
    %print progress and status without flooding the command window
    while nbytes > 0
        fprintf('\b')
        nbytes = nbytes - 1;
    end
    nbytes = fprintf('processing %.3f of %d, QP status: %s', Ts*k, Tmax, yalmiperror(err));
    
    %update (measured) output
    y(:,k) = C*x(:,k) + w(:,k);

    if(k >= Tini+1)
        %update controller
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
        %to prevent infeasible solutions initially, let the system run in
        %open loop until enough data is gathered to construct U_ini, Y_ini
       du(:,k) = 1*rand(nu, 1); %u(:,k) = 1*rand(nu, 1);
    end

    %update system dynamics
    u(:,k) = u(:,max(1,k-1)) + du(:,k);
    
    %during the initial part of the simulation when random inputs are
    %generated, the following operations limit the u to the interval [0,1]
    %where the duty-cycles are constrained
    
    %after k>=Tini+1 these operations will have no effect, because the
    %offset-free SPC solvers enforces the constraints
    
    if u(1,k) >=1
        u(1,k)=1;
    end
    if u(1,k) <=0
        u(1,k)=0;
    end
    if u(2,k) >=1
        u(2,k)=1;
    end
    if u(2,k) <=0
        u(2,k)=0;
    end
    
    x(:,k+1) = A*x(:,k) + B*u(:,k) + [0;0;0;0;1]*d(:,k);
    
    
          
end

%% Display simulation results

figure();
ax1 = subplot(311);
stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference');
hold on;
stairs(t, y(1,:), 'r', 'DisplayName', 'iSPC');
axis([0 t(end) 0 5.5]);

ylabel('Output current [A]');
legend;

ax2 = subplot(312);
stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference');
hold on;
stairs(t, y(1,:), 'r', 'DisplayName', 'iSPC');

ylabel('Output current [A]');
axis([0 t(end) 4.5 5.5]);

ax3 = subplot(313);
hold on;
stairs(t, u(1,:), 'r');
stairs(t, u_ss(1)*ones(1,length(t)), '--k');
stairs(t, u(2,:), 'b');
stairs(t, u_ss(2)*ones(1,length(t)), '--k');

ylabel('Duty-cycles ');





