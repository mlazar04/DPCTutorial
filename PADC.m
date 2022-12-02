function [sys, ymin, ymax, umin, umax, Mx, Nx, Mu, Nu] = PADC()

%% Author: Xu Duo

%% Amplifier circuit parameters
Vbus = 360;
Rm = 10;
Lm = 20e-3;
R = 60.2e-6;
L = 44e-6;
C = 0.4e-6;

%% Averaged continuous-time state-space model
Ac = [-R/L -1/L 0 0 +R/L;
      1/C 0 0 0 -1/C;
      0 0 -R/L -1/L -R/L;
      0 0 1/C 0 1/C;
      R/Lm 1/Lm -R/Lm -1/Lm -(2*R+Rm)/Lm];
Bc = [Vbus/L 0;
      0 0;
      0 Vbus/L;
      0 0;
      0 0];
Cc = [0 0 0 0 1]; 

sys = ss(Ac, Bc, Cc, 0);

%% Output current constraints
ymin = -30;
ymax = 30;

%% Duty-cycle constraints
umin = [0;0];
umax = [1;1];

%% Variables for computing steady-state values
Mx = [1;Rm/2;-1;-Rm/2;1];
Nx = [0;Vbus/2;0;Vbus/2;0];
Mu = [Rm/2/Vbus;-Rm/2/Vbus];
Nu = [1/2;1/2];
end

