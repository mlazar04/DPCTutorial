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

%% computation of steady-state equilibrium values
x_ss = Mx*5+Nx;
u_ss = Mu*5+Nu;

Wvar = 0.1; %noise variance for offline measured data
wvar = 0.1; %noise variance for online measured data


%% Controller parameters
N = 15;               % Prediction horizon
Q = 10;               % Output cost weight
R = 0.1*eye(2);       % Input cost weight
Tini = 12;            % Window of past input data used online for initializing the predicted trajectories




%% Generate data
T = 200;  %the larger the data length, the better the robust performance; complexity of DeePC increases; 
U = 0.8*idinput([T+N+Tini, nu], 'PRBS', [0, 1], [0, 1])';
l1=1*1e+5;                % Weight of g regularization (zero if there is no noise)
l2=1*1e+5;                %lambda 2 (zero if there is no noise or disturbance)
l3=1*1e+5;                %lambda 3 (zero if there is no noise or disturbance)

X = zeros(n, size(U,2));
Y = zeros(ny,size(U,2));

%sensor noise
W = Wvar*randn(ny,size(U,2));

for k=1:size(U,2)
    Y(:,k) = C*X(:,k) + W(:,k);
    if(k < size(U,2))
        X(:,k+1) = A*X(:,k) + B*U(:,k);
    end
end

disp(['SNR output noise: ', num2str(snr(Y, W))]);

%% Define DeePC Controller
Up = zeros(Tini*nu, T);
Uf = zeros(N*nu, T);
Yp = zeros(Tini*ny, T);
Yf = zeros(N*ny, T);

for i = 1:Tini
    Up((i-1)*nu+1:i*nu, :) = U(:, i  :i+T-1);
    Yp((i-1)*ny+1:i*ny, :) = Y(:, i+1:i+T);
end

for i = 1:N
    Uf((i-1)*nu+1:i*nu, :) = U(:, i+Tini  :i+Tini+T-1);
    Yf((i-1)*ny+1:i*ny, :) = Y(:, i+Tini+1:i+Tini+T);
end

%computing matrix weight for least squares regularization of g
PI = pinv([Up; Yp; Uf])*[Up; Yp; Uf];

Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);

%Build YALMIP DeePC problem
Uss=[];
for i=1:1:N
    Uss=[Uss; u_ss];
end

u = sdpvar(nu*N,1);
y = sdpvar(ny*N,1);
ref = sdpvar(ny*N,1);
g = sdpvar(T, 1); 
sigma_u = sdpvar(nu*Tini,1);
sigma_y = sdpvar(ny*Tini,1);
u_ini = sdpvar(Tini*nu, 1);
y_ini = sdpvar(Tini*ny, 1);

% Define objective function
objective = (y-ref)'*Omega*(y-ref)+(u-Uss)'*Psi*(u-Uss);  % tracking DeePC cost
objective = objective + l1*(g'*(eye(length(PI))-PI)'*(eye(length(PI))-PI)*g) + l2*(sigma_u'*sigma_u) + l3*(sigma_y'*sigma_y);  %regularization terms

%above the term g'*(eye(length(PI))-PI)'*(eye(length(PI))-PI)*g can be
%replaced by g'*g to obtain the original DeePC regularization, which is
%however not consistent

% Build constraints
constraints = [u_ini==Up*g+1*sigma_u, y_ini==Yp*g+1*sigma_y, y==Yf*g, u==Uf*g];  % DeePC equality constraints

%if there is no disturbance / noise multiply sigma_u and sigma_y by zero to
%remove the softening of the equality constraints

for k = 1:N   % and the bounds on the inputs and outputs
    constraints = [constraints, ymin<=y(ny*(k-1)+1:ny*k)<=ymax, umin<=u(nu*(k-1)+1:nu*k)<=umax];
end

Parameters = {u_ini, y_ini, ref};
Outputs = {u, y};

%% Quadprog is a standard Matlab QP solver; Mosek is numerically more consistent; other solvers can be used with YALMIP
options = sdpsettings('solver', 'quadprog', 'verbose', 0, 'debug', 0);
controller = optimizer(constraints, objective, options, Parameters, Outputs);

%% Initialize Simulation
Tmax = 600*Ts;
t = 0:Ts:Tmax;
simLen = size(t,2);
x0 = zeros(5,1); %initial condition

%reference and disturbance sequences. Note that reference must be N samples
%longer to accommodate for the last predictions
r = zeros(ny, simLen+N);
r(:,5e-5/Ts:end) = 5;
d = zeros(ny, simLen);
d(:,6e-4/Ts:1e-3/Ts) = 1*-0.025;%active (multiply by 1) or inactive disturbance (multiply by 0)



y = zeros(ny, simLen);
u = zeros(nu, simLen);
x = zeros(n, simLen+1);
x(:,1) = x0;

w = wvar*randn(ny, simLen);

%% Simulation
nbytes = fprintf('time: 0 of %d', Tmax);
err = 0;

for k = 1:simLen
    
    %Print current time and status without flooding the command window
    while nbytes > 0
        fprintf('\b')
        nbytes = nbytes - 1;
    end
    nbytes = fprintf('processing %.3f of %d, QP status: %s', Ts*k, Tmax, yalmiperror(err));
    
    %Compute new (measured) output
    y(:,k) = C*x(:,k) + w(:,k);

    if(k >= Tini+1)
        %Update control law
        Rk = r(:, k+1:k+N);
        Rk = Rk(:);
        U_ini = u(:, k-Tini:k-1);
        U_ini = U_ini(:); %flatten vector
        Y_ini = y(:,k-Tini+1:k);
        Y_ini = Y_ini(:);
        
        [Sol, err] = controller({U_ini, Y_ini, Rk});
        Uk = Sol{1};
        Yk = Sol{2}; %Yk can be used for evaluation, is not stored
        
        u(:,k) = Uk(1:nu);
    else
        % To prevent initial infeasibility, let the system run in open loop
        % until enough data is gathered to build U_ini, Y_ini
        u(:,k) = 1*rand(nu, 1);
    end

    % Update system
    x(:,k+1) = A*x(:,k) + B*u(:,k) +[0;0;0;0;1]*[d(:,k)];
end

%% Display simulation results

figure();
ax1 = subplot(311);
stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference');
hold on;
stairs(t, y(1,:), 'r', 'DisplayName', 'DeePC');

ylabel('Output current [A]');
legend;
axis([0 t(end) 0 5.5]);
ax2 = subplot(312);
stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference');
hold on;
stairs(t, y(1,:), 'r', 'DisplayName', 'DeePC');

ylabel('Output current [A]');
axis([0 t(end) 4.98 5.02]);

ax3 = subplot(313);
hold on;
stairs(t, u(1,:), 'r');
stairs(t, u_ss(1)*ones(1,length(t)), '--k');
stairs(t, u(2,:), 'b');
stairs(t, u_ss(2)*ones(1,length(t)), '--k');

ylabel('Duty-cycles');
xlabel('Time [s]');
