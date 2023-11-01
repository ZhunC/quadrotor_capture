clc
clear


syms x y z vx vy vz phi theta psi w1 w2 w3 u1 u2 u3 u4 m l sig real
syms Ix Iy Iz

% definitions
I = diag(sym([Ix Iy Iz]));
pos = [x, y, z]';
vel = [vx, vy, vz]';
a = [phi, theta, psi]';
w = [w1, w2, w3]';
states = [pos; a; vel; w];
u = [u1, u2, u3, u4]';
g = 9.81;

% rotation matrices
Tinv = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
    0, cos(phi), -sin(phi);...
    0, sin(phi)/cos(theta), cos(phi)/cos(theta)];

R3 = [cos(psi), -sin(psi), 0;...
    sin(psi), cos(psi), 0;...
    0, 0, 1];
R2 = [cos(theta), 0, sin(theta);...
    0, 1, 0;...
    -sin(theta), 0, cos(theta)];
R1 = [1, 0, 0;...
    0, cos(phi), -sin(phi);...
    0, sin(phi), cos(phi)];
Rce = R3*R2*R1;
Rec = transpose(Rce);


%% state space equations
f1 = vel;
f2 = Tinv * w;
f3 = Rce*(u1+u2+u3+u4)/m*[0;0;1]+[0;0;-g];  % assume no external force
f4 = inv(I)*[(u2-u4)*l, (u3-u1)*l, (u1-u2+u3-u4)*sig]'-cross(w,I*w); % no ext torque

% map body axis to space axis in f4
%% taking jacobians
f = [f1; f2; f3; f4];
A = jacobian(f, states);
B = jacobian(f, u);


%% evaluating at i.c.
syms x1 y1 z1 real
states_0 = [x1;y1;z1;0;0;0;0;0;0;0;0;0];
u0 = [m*g/4;m*g/4;m*g/4;m*g/4];
Anum = subs(A,[states;u],[states_0;u0]);
Bnum = subs(B,[states;u],[states_0;u0]);
%% check controllability
Anum_double = double(Anum);
% assume values for the constants
% m = 0.5; %kg
% l = 0.1; %m
% Iz = 0.02;
% Iy = 0.02;
% Iz = 0.05;
% sig = 0.9
Bnum_double = double(subs(Bnum,[m;l;Ix;Iy;Iz;sig],[0.5; 0.1; 0.02; 0.02; 0.05; 0.9]));

Co = ctrb(Anum_double,Bnum_double);
Co_rank = rank(Co);

%% open loop poles

C = [eye(3),zeros(3,9)];
D = zeros([3,4]);
sys = ss(Anum_double, Bnum_double, C, D);
poles = pole(sys);


%%
% assume I want all poles at -2
% p = -2*ones([12,1]);
p = -13:-2;

K = place(Anum_double, Bnum_double, p);

% the above thing is working

%% verify

figure(1)
step(sys)
hold on
Acl = Anum_double - Bnum_double*K;
sys2 = ss(Acl, Bnum_double, C, D);
figure(1)
step(sys2)

xlim([0 10])




