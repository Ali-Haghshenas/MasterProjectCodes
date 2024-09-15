%%%%%%%%% Parameters of OneCell %%%%%%%%%
dt = 0.05; %[sec.]
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

C = [0 -1 -1];

D = -R0;

sys_c = ss(A,B,C,D);

%%%%%%%%% Discrete time System %%%%%%%%%%
%%%%%%%%% OneCell equation:

%Initial condition (OneCell)
X0 = [0.6;0;0];

sys_d = c2d(sys_c,dt,'zoh');
A_d = sys_d.A;
B_d = sys_d.B;

%%%%%%%%%% Pack equation:
%%%%%%%%% Parameters of Battery Pack %%%%%%%%%
taw_1_pack = 808*[1;1;1;1]; %[808;800;802;805]; %[sec.]
taw_2_pack = 36*[1;1;1;1]; %[36;37;38;35]; %[sec.]
C_max_pack = 8640*[1;1;1;1]; %[8640;8641;8639;8642]; %[As]==[Farad] eqals to 2.4269 [Ah]
R_D1_pack = 199e-3*[1;1;1;1]; %[ohm]
R_D2_pack = 25e-3*[1;1;1;1]; %[ohm]
%We use 4 Li-ion battery cells connected serially.
R0_pack = [R0;R0;R0;R0];

A1 = [0 0 0;
    0 1/taw_1_pack(1,1) 0;
    0 0 1/taw_2_pack(1,1)];

B1 = [-eta/C_max_pack(1,1);
    1/(taw_1_pack(1,1)/R_D1_pack(1,1));
    1/(taw_2_pack(1,1)/R_D2_pack(1,1))];

A2 = [0 0 0;
    0 1/taw_1_pack(2,1) 0;
    0 0 1/taw_2_pack(2,1)];

B2 = [-eta/C_max_pack(2,1);
    1/(taw_1_pack(2,1)/R_D1_pack(2,1));
    1/(taw_2_pack(2,1)/R_D2_pack(2,1))];

A3 = [0 0 0;
    0 1/taw_1_pack(3,1) 0;
    0 0 1/taw_2_pack(3,1)];

B3 = [-eta/C_max_pack(3,1);
    1/(taw_1_pack(3,1)/R_D1_pack(3,1));
    1/(taw_2_pack(3,1)/R_D2_pack(3,1))];

A4 = [0 0 0;
    0 1/taw_1_pack(4,1) 0;
    0 0 1/taw_2_pack(4,1)];

B4 = [-eta/C_max_pack(4,1);
    1/(taw_1_pack(4,1)/R_D1_pack(4,1));
    1/(taw_2_pack(4,1)/R_D2_pack(4,1))];

A_pack = [A1 zeros(3,9);zeros(3,3) A2 zeros(3,6);zeros(3,6) A3 zeros(3,3);zeros(3,9) A4];
B_pack = [B1;B2;B3;B4];
C_pack = [C C C C];
D_pack = -R0_pack(1,1)-R0_pack(2,1)-R0_pack(3,1)-R0_pack(4,1);

sys_c_pack = ss(A_pack,B_pack,C_pack,D_pack);
sys_d_pack = c2d(sys_c_pack,dt,'zoh');
A_pack_d = sys_d_pack.A;
B_pack_d = sys_d_pack.B;
% A_pack = kron(eye(4),A_d);
% B_pack = [B_d;B_d;B_d;B_d];

%Initial Condition (BatteryPack)
X10 = [0.6;0;0];
X20 = [0.8;0;0];
X30 = [0.7;0;0];
X40 = [0.9;0;0];
Zeta0 = [X10;X20;X30;X40];