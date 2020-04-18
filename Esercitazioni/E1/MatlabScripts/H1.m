% Computational Thermofluid dynamics
% Homework No.1 - Determination of the heat flow through a fin.

% Author: Davide Armenante
% See the quantities in the attacched relation.

% k: thermal conductivity
% t: fin thickness
% L: fin length
% h: convection heat transfer coefficient

clear all;
close all;
clc;

% Variables definition
% Geometric propertyes
t       = 20;      % mm
L       = 200;      % mm

% Material thermal propertyes
k       = 50;       % W/(m K)
h       = 500;      % W/(m^2 K)

% Temperature definition
Tinf    = 25;       % Celsius Degree
Tb      = 200;      % Celsius Degree

% Number of Finite Volumes
N       = 100;
%Conversion
t   = t * 10^-3;
L   = L * 10^-3;
Ac  = t;            % Fin cross-section of unitary depth

%% Analytic solution
m           = sqrt((2*h)/(k*t));
qfExact     = sqrt(2*h*k*t)*(sinh(m*L)+(h/(m*k))*cosh(m*L))/...
                (cosh(m*L)+(h/(m*k))*sinh(m*L))*(Tb - Tinf);

%% Array preallocation
AE   = zeros(N,1);  % East coefficient
AW   = zeros(N,1);  % West coefficient
AP   = zeros(N,1);  % Central coefficient
SP   = zeros(N,1);  % Source
T    = zeros(N,1);  % Temperature

%% 
%Fin distretization
Dx        = L/N;                       % Length of each FV
Xf        = linspace(0,L,N+1);         % Position of CV faces
inode     = 1:N;                       % Index of all nodes
Xc(inode) = (Xf(inode)+Xf(inode+1))/2; % Position of CV centroids
DAs       = 2*Dx;                      % External surface for each CV

AE(1:N-1) = (k/Dx)*Ac;
AW(2:N)   = (k/Dx)*Ac;
AP(:)     = -(AW+AE)-DAs*h;
SP(:)     = - DAs*h*Tinf;

% Boundary condition
% Fixed Temperature Tb at the base
AP(1)     = AP(1)-(2*k/Dx)*Ac;
SP(1)     = SP(1)-(2*k/Dx)*Ac*Tb;

% Convection at the tip
AP(N)     = AP(N)-h*Ac;
SP(N)     = SP(N)-h*Tinf*Ac;

%% Sparse diagonal matrix assembly
A   = spdiags([circshift(AW,[-1 0]) AP circshift(AE,[1 0])],...
        [-1 0 1], N, N);

% Solution
T   = A\SP;

%% Conductive flux an relative error
qfNum   = -2*k*Ac*(T(1)-Tb)/Dx;
diff    = (abs(qfNum-qfExact))/qfExact*100;

disp(['Theoretical heat flux = ', num2str(qfExact),' [W/m]']);
disp(['Difference Num.-Analyt. = ', num2str(diff),' %']);

%% Plot
plot(Xc,T,'-o','MarkerFaceColor','b','MarkerSize',5);
xlim([0, Xc(end)*1.1])
ylim([T(end)*.9, T(1)*1.1])
xlabel('x [m]');
ylabel('T [C]');
title('Temperature profile along the fin');

%% Save data
Xc          = Xc';
TempData    = table(Xc,T);
save('Data/qf_Exact1.mat','qfExact');
save('Data/aratio1.mat','TempData');