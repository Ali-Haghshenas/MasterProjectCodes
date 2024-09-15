%% Main Section %%
%%%%%%%%% Parameters of OneCell %%%%%%%%%
T_s = 0.001; %0.000008 Sampling Time [sec.]
taw_1 = 808; %[sec.]
taw_2 = 36; %[sec.]
eta = 1; %We consider that the battery is charged and discharged with efficiency 1
C_max = 8640; %[As]==[Farad] eqals to 2.4269 [Ah]
R_D1 = 199e-3; %[ohm]
R_D2 = 25e-3; %[ohm]
R0 = 43e-3; %[ohm]

%%%%%%%%% Continuous time System %%%%%%%%%%
A = [0 0 0;
    0 1/taw_1 0;
    0 0 1/taw_2];

B = [-eta/C_max;
    1/(taw_1/R_D1);
    1/(taw_2/R_D2)];

C = [0.5672 -1 -1];

D = -R0;

bias = 3.551;

sys_c = ss(A,B,C,D);

%%%%%%%%% Discrete time System %%%%%%%%%%
%%%%%%%%% OneCell equation:
sys_d = c2d(sys_c,T_s,'zoh');
A_d = sys_d.A;
B_d = sys_d.B;
C_d = sys_d.C;
D_d = sys_d.D;

%Initial condition (OneCell)
% X0 = [0.6;0;0];

%%%%%%%%% Observer:
rank_obsv = rank(obsv(A_d,C_d));
L = acker(A_d',C_d',[0.001,0.001,0.001])';

%% MPC Section %%
%%%%%%%%% Augmented System %%%%%%%%%%
A_aug = [A_d zeros(3,1);C_d*A_d 1];
B_aug = [B_d;C_d*B_d+D_d];
C_aug = [zeros(1,3) 1];
D_aug = D_d;

%%%%%%%%% MPC Simulation %%%%%%%%%%
%%%%%%%%% Controller Parameters:
N_p = 1000;        % Prediction horizon  
N_c = 1000;         % Control horizon 
R_s_bar = ones(N_p,1);
R_bar = eye(N_c);
I_bar = eye(N_c);

%%%%%%%%% Augmented Observer:
L_aug = acker(A_aug', C_aug',[0.001,0.001,0.001,0.001])';

%%%%%%%%% MPC gains:
[Phi,F] = mpcgain(A_aug,B_aug,C_aug,N_c,N_p);
 H = (Phi'*Phi+R_bar);
 Phi_R = Phi'*R_s_bar;
 Phi_F = Phi'*F;
 K_r = [1 zeros(1,N_c-1)]*inv(H)*Phi_R;
 K_mpc = [1 zeros(1,N_c-1)]*inv(H)*Phi_F;

%%%%%%%%% Simulation Parameters:
sim_time = 1;   % Total simulation time [seconds]  
num_steps = round(sim_time / T_s);
Delta_u_history = zeros(num_steps, 1); % For storing the genearted input with mpc  
Y_history = zeros(num_steps, 1); % For storing the output  

%%%%%%%%% State Vector Initialization  
X0 = [0.6; 0; 0];  % Initial condition  
X_k = X0;          % Initial state
x0 = [0;0;0;0];    % Initial state of Augmented system
u_k = 0;           % Initial input
r_k = 3.85;        % Reference output signal (initial)
x_hat_k = zeros(size(A_aug,1),1);  % State Estimate (initial)
Delta_u_k = 0;
Y_k = C_d * X_k + D_d * u_k + bias; %If the plant and the model are not the same, please change this to real values of plant

%%%%%%%%% Reference Output Signal:
R_s = R_s_bar*r_k;         % Reference Output (Y_k that reported to controller from real plant)

%%%%%%%%% Constraints:
u_max = 4.8;        % Maximum control input  
u_min = -2.4;       % Minimum control input  
Y_max = 4.2;        % Maximum of battery output
Y_min = 0;       % Minimum of battery output
y_min = Y_min * ones(N_p, 1);
y_max = Y_max * ones(N_p, 1);
Delta_u_max = 0.1;
Delta_u_min = -0.1;
Delta_U_max = ones(N_c,1)*Delta_u_max;
Delta_U_min = ones(N_c,1)*Delta_u_min;
W_mmax = -D_aug*Delta_u_min*ones(N_p,1); % -D_aug*Delta_u_max*onesones(N_p,1)<=-W<=-D_aug*Delta_u_min*onesones(N_p,1)
W_pmax = D_aug*Delta_u_max*ones(N_p,1); % D_aug*Delta_u_min*onesones(N_p,1)<=W<=D_aug*Delta_u_max*onesones(N_p,1)

for k_i = 0:num_steps-1
    % Reference signal u(k)
    r_k1 = 3.85; %cte signal
    % Input signal u(k)
    u_k1 = 2.4; %cte signal
    Delta_u_k1 = u_k1 - u_k;

    % Update states
    X_k1 = A_d * X_k + B_d * u_k;
    Delta_X_k1 = X_k1 - X_k;

    Y_k1 = C_d * X_k1 + D_d * u_k1 + bias;

    x_k1 = [Delta_X_k1' Y_k1-D_d*Delta_u_k1]';

    x_hat_k1 = A_aug * x_hat_k + B_aug * Delta_u_k + L_aug * (Y_k - C_aug * x_hat_k - D_aug * Delta_u_k);

    % Contraints Matrices
    N_1 = [u_max-u_k;-u_min+u_k];
    C_1 = [1 zeros(1,N_c-1)];
    M_1 = [C_1;-C_1];
    N_2 = [Delta_U_max;-Delta_U_min];
    M_2 = [I_bar;-I_bar];
    N_3 = [y_max-F*x_hat_k1+W_mmax;-y_min+F*x_hat_k1+W_pmax];
    M_3 = [Phi;-Phi];
    M = [M_1;M_2;M_3];
    gamma = [N_1;N_2;N_3];
    
%     H = (Phi'*Phi+R_bar);
%     Phi_R = Phi'*R_s_bar;
%     Phi_F = Phi'*F;
    
    % Constraints gain
%     lambda = 2 * inv(M*inv(H)*M') * (M * inv(H) * (Phi_R * r_k1 - Phi_F * x_hat_k1) - gamma); %Just for constrainted problem

    % Delta U mpc
%     Delta_U = inv(H) * (Phi_R * r_k1 - Phi_F * x_hat_k1 - M'*lambda/2); % By implementing constraints
    Delta_U = inv(H) * (Phi_R * r_k1 - Phi_F * x_hat_k1); % Without any constraints

    % Receding Horizon
    Delta_u_k1_mpc = [1 zeros(1,N_c-1)] * Delta_U;

    Delta_u_history(k_i+1,1) = Delta_u_k1_mpc;
    Delta_u_k1 = Delta_u_k1_mpc;

    % Update previous variables for next sample
    u_k = u_k1;
    X_k = X_k1;
    Y_k = Y_k1;
    x_hat_k = x_hat_k1;
    Delta_u_k = Delta_u_k1;
end

k = 1:1:num_steps;
figure;  
stairs(k , Delta_u_history');
xlabel('Time (samples)');  
ylabel('\Delta u(k_i)');  
title('\Delta u(k_i) Results');  
grid on;