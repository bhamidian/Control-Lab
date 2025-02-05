clear
close all
clc
format compact

%% Cart-Pendulum Linearization using Symbolic Math Toolbox
syms F real
syms x [4, 1] real
mc = 0.5;
mp = 2;
L = 0.2;
g = 9.81; 

% Define nonlinear functions
f1 = ((mc + mp) * g * sin(x3) - mp * L * x4^2 * sin(x3) * cos(x3)) / ((mc + mp * sin(x3)^2) * L);
g1 = -cos(x3) / ((mc + mp * sin(x3)^2) * L);
f2 = (-mp * g * sin(x3) * cos(x3) + mp * L * x4^2 * sin(x3)) / (mc + mp * sin(x3)^2);
g2 = 1 / (mc + mp * sin(x3)^2);

% Define the nonlinear state equations
f = [
    x2;
    f2 + g2 * F;
    x4;
    f1 + g1 * F
];

% Linearization
x_eq = [0; 0; pi; 0]; % Equilibrium state
F_eq = 0;            % Equilibrium input

% Compute Jacobians
A = jacobian(f, x);
B = jacobian(f, F);

% Substitute equilibrium points
A = subs(A, [x; F], [x_eq; F_eq]);
B = subs(B, [x; F], [x_eq; F_eq]);

% Output C and D matrices for linearized system
C = [1, 0, 0, 0];
D = 0;

% Simplify results
A = double(simplify(A));
B = double(simplify(B));

% Display results
disp('Linearized System Matrices:');
disp('A ='); disp(A);
disp('B ='); disp(B);
disp('C ='); disp(C);
disp('D ='); disp(D);

Sys_TF = tf(ss(A, B, C, D))
eig(A)
