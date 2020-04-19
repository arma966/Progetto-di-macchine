%% Plot different Temperature profile for different fin thickness

clear all
load('Data/aratio1.mat')
load('Data/qf_Exact1.mat')

plot(TempData.Xc,TempData.T,'-o','MarkerFaceColor','b','MarkerSize',5);
xlabel('x [m]');
ylabel('T [°C]');
title('Temperature profile along the fin');

hold on

text(TempData.Xc(end/2)*.8,TempData.T(end/2)*1.4,strcat("qf = ",num2str(round(qfExact,0))," [W/m]"))
clear TempData

load('Data/aratio2.mat')
load('Data/qf_Exact2.mat')
plot(TempData.Xc,TempData.T,'-v','MarkerFaceColor','r','MarkerSize',5);
text(TempData.Xc(end/2)*.8,TempData.T(end/2)*1.4,strcat("qf = ",num2str(round(qfExact,0))," [W/m]"))
legend('L/t = 10','L/t = 2')

hold off

%% Save
print('../../LaTex/figH/TempProf','-deps')