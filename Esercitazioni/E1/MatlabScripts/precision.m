% Computational Thermofluid dynamics
% Homework No.1 - Determination of the heat flow through a fin.

% Determination of the relative error as a function of the number of finite
% elements.

% Author: Davide Armenante
% See the quantities in the attacched relation.

% Definitions
% k: thermal conductivity
% t: fin thickness
% L: fin length
% h: convection heat transfer coefficient

clear all;
close all;
clc;

% Variables definition
% Geometric propertyes
t       = 20;           % mm
L       = 200;          % mm

% Material thermal propertyes
k       = 50;           % W/(m K)
h       = 500;          % W/(m^2 K)

% Temperature definwition
Tinf    = 25;           % Celsius Degree
Tb      = 200;          % Celsius Degree

%Conversion
t       = t * 10^-3;
L       = L * 10^-3;

Ac      = t;            % Fin cross-section of unitary depth

%% Analytic solution
m           = sqrt((2*h)/(k*t));
qfExact     = sqrt(2*h*k*t)*(sinh(m*L)+(h/(m*k))*cosh(m*L))/...
                (cosh(m*L)+(h/(m*k))*sinh(m*L))*(Tb - Tinf);
%%
for i = 1:10
    disp(['Passo: ', int2str(i)]);
    
    N    = 10*i         %Number of fin subdivision
    
    % Array preallocation
    AE   = zeros(N,1);  % East coefficient
    AW   = zeros(N,1);  % West coefficient
    AP   = zeros(N,1);  % Central coefficient
    SP   = zeros(N,1);  % Source
    T    = zeros(N,1);  % Temperature

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

    % Sparse diagonal matrix assembly
    A   = spdiags([circshift(AW,[-1 0]) AP circshift(AE,[1 0])],...
            [-1 0 1], N, N);

    % Solution
    T = A\SP;

    % Conductive flux
    qfNum       = -2*k*Ac*(T(1)-Tb)/Dx;
    diff        = (abs(qfNum-qfExact))/qfExact*100;
    Error(i)    = diff;
    
    disp(['Difference Num.-Analyt. = ', num2str(diff),' %']);
end

%% Plot
xaxis = 10:10:N;
plot(xaxis,Error, '-o','MarkerFaceColor','b','MarkerSize',5)
xlim([0, 1.1*N])
ylim([-.2, Error(1)*1.05])
ylabel({'$\epsilon \%$'},'Interpreter','latex');
xlabel('N');
title('Relative error')
% Uncomment to insert the plot in the final issue
print('../../LaTex/figH/Error','-deps')