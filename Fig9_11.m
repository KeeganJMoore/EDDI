%% 
clc,clear,close all
tic
load("AllData_Processed_SN.mat")

%%
clc,close all
[(1:length(Data.MaxForces))' Data.MaxForces]

mass = 0.815;
Expi = 18; % Use 18 for ID, then 12, 13, 19, 25 for generalization

time = Data.Time;
Acc = Data.Acc(:,Expi);
Vel = Data.Vel(:,Expi);
Disp = Data.Disp(:,Expi);
Force = Data.Force(:,Expi);
MaxF = Data.MaxForces(Expi);
F1 = 2; F2 = 20;
Disp = MoDAL.IFD(time,Disp,F1,F2);
Vel = MoDAL.IFD(time,Vel,F1,F2);
Acc = MoDAL.IFD(time,Acc,F1,F2);

MoDAL.PlotForce(time,Force,"fontSize",12);
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

MoDAL.PlotTSWTFT(time,1000*Disp,0,20,"label",'Dispmm','title', ...
    ['Displacement Response, Force = ' num2str(round(MaxF)) ' N'], ...
    "fontSize",14,"timeEnd",45)
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%%

% Find Time of Peak Velocity
[maxV,vMaxInd] = max(abs(Vel));
timeID = time(vMaxInd:end);
VelID = Vel(vMaxInd:end);
DispID = Disp(vMaxInd:end);

[~,idx] = findpeaks(-abs(DispID));


KEId = 1/2*mass*VelID.^2;
R = KEId(1)-KEId(idx);

% Plot Kinetic Energy and KE when displacement is zero.
figure
plot(timeID,KEId,'k'); hold on; plot(timeID(idx),KEId(idx),'go','MarkerSize',5);
hold on
tin = 4; tfin = 6; yin = 0; yfin = 0.008;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
xlim([0 45])
xlabel('Time [s]'); ylabel('Energy [J]');
legend({'Experimental KE, T(t)', 'T(t=\gamma_i)'},Location='southeast'); % , 'Interpreter', 'latex',Location='southeast'
set(gca,'fontsize',14)

p2 = axes;
p2.Position = [0.40 0.35 0.45 0.5];
plot(timeID,KEId,'k'); hold on; plot(timeID(idx),KEId(idx),'go','MarkerSize',5);
xlim([tin tfin])
annotation('Arrow','Position',[0.235,0.2663,0.0986,0.0620],'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%%

% Compute Damping Coefficients
Pow = 1;
A = zeros(length(idx),Pow); 
for m = 2 : Pow+1
    xdot = VelID.^ m;
    At = cumtrapz(timeID,xdot);
    A(:,m-1) = At(idx);
end
x_xdot = DispID.^ 2.*VelID.^2;
RR = cumtrapz(timeID,x_xdot);
A1 = RR(idx);

xdot2_xdot = VelID.^4;
RR2 = cumtrapz(timeID,xdot2_xdot);
A2 = RR2(idx);

Am1 = [A A1 A2];
% Am2 = [A A2];

dampModel1 = Am1\R

% Compute Model Forces & Energies
dissForce1 = dampModel1(1)*VelID+dampModel1(2)*DispID.^2.*VelID+dampModel1(3)*VelID.^3;

dissEnergy1 = cumtrapz(timeID,VelID.*dissForce1);

MEEst_Model1 = KEId(1)-dissEnergy1;

% Plot dissipated energies
figure
plot(timeID(idx),R,'ko',timeID,dissEnergy1,'r')
hold on
plot([4 6],[0.02 0.02],'color',[0.7 0.7 0.7])
plot([4 6],[0.026 0.026],'color',[0.7 0.7 0.7])
plot([4 4],[0.02 0.026],'color',[0.7 0.7 0.7])
plot([6 6],[0.02 0.026],'color',[0.7 0.7 0.7])
xlim([0 45])
xlabel('Time [s]')
ylabel('Energy [J]')
% title("Comparison of Dissipated Energies")
legend('Experimental Dissipated Energy','Estimated Dissipated Energy','location', ...
    'southeast')

p2 = axes;
p2.Position = [0.4375 0.39048 0.45 0.42];
plot(timeID(idx),R,'ko',timeID,dissEnergy1,'r')
xlim([4 6])
annotation('Arrow','Position',[0.234 0.688 0.092857 -0.02619],'color',[0.7 0.7 0.7]) 
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

% Plot mechanical energies
figure
plot(timeID,KEId,'k',timeID,MEEst_Model1,'r')
hold on
plot([4 6],[0 0],'color',[0.7 0.7 0.7])
plot([4 6],[0.01 0.01],'color',[0.7 0.7 0.7])
plot([4 4],[0 0.01],'color',[0.7 0.7 0.7])
plot([6 6],[0 0.01],'color',[0.7 0.7 0.7])
ylim([0 0.03]); xlim([0 45])
xlabel('Time [s]');ylabel('Energy [J]')
% title("Identified Mechanical Energy")
% title(['Force = ' num2str(MaxF) ' N'])
legend('Experimental KE','Estimated ME', 'Interpreter', 'latex','position', ...
    [0.50009 0.13801 0.34732 0.15595])
set(gca,'fontsize',12)

p2 = axes;
p2.Position = [0.40 0.375 0.45 0.5];
plot(timeID,KEId,'k',timeID,MEEst_Model1,'r')
xlim([4 6])
annotation('Arrow','Position',[0.235 0.32619 0.10357 0.059524], ...
    'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%%
% Fstart = 30000;
Fskip = 10;
potEnergy1 = MEEst_Model1 - KEId;

% Plot potential energy
figure
plot(timeID,potEnergy1,'k')
hold on
plot([4 6],[-0.002 -0.002],'color',[0.7 0.7 0.7])
plot([4 6],[0.01 0.01],'color',[0.7 0.7 0.7])
plot([4 4],[-0.002 0.01],'color',[0.7 0.7 0.7])
plot([6 6],[-0.002 0.01],'color',[0.7 0.7 0.7])
legend('Estimated Potential Energy','position',[0.41429 0.1369 0.44643 0.045238]);
xlim([0 45])
xlabel('Time [s]'); ylabel('Energy [J]'); set(gca,'fontsize',12)

p2 = axes;
p2.Position = [0.425 0.375 0.45 0.5];
plot(timeID,potEnergy1,'k')
xlim([4 6])
h = annotation('Arrow','Position',[0.235 0.39762 0.13464 0.028572], ...
    'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

% consForce1 = diff(potEnergy1)./diff(DispID);

L = 2*KEId-MEEst_Model1;
F1Est = gradient(L)./gradient(DispID) - mass*Acc(vMaxInd:end);      
consForce1 = -F1Est;
consForce1 = detrend(smoothdata(consForce1,'movmedian',200));     % Modified by CL

figure
p1 = axes;
plot(timeID(1:Fskip:end-1),consForce1(1:Fskip:end),'k')
hold on
plot([4 6],[-7 -7],'color',[0.7 0.7 0.7])
plot([4 6],[7 7],'color',[0.7 0.7 0.7])
plot([4 4],[-7 7],'color',[0.7 0.7 0.7])
plot([6 6],[-7 7],'color',[0.7 0.7 0.7])
xlabel('Time [s]'); ylabel('Force [N]'); xlim([0 45]); ylim([-12 12])
legend("Estimated $K$", 'Interpreter', 'latex','Location','southeast')
set(gca,'fontsize',12)

p2 = axes;
p2.Position = [0.40 0.6 0.45 0.25];
plot(timeID(1:Fskip:end-1),consForce1(1:Fskip:end),'k')
set(gca,'fontsize',14)
xlim([4 6]);ylim([-7 7])
annotation('Arrow','Position',[0.235,0.6520,0.1086,0.00], ...
    'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%% Identify model for conservative force
dID = DispID;

x = dID;
theta = [x  x.^2  x.^3  x.^4  x.^5];
tic
K = theta\consForce1

%%
kModel = [K(1); K(2); K(3); K(4); K(5)]

ModFc = kModel(1)*dID + kModel(2)*dID.^2 +  kModel(3)*dID.^3+ kModel(4)*dID.^4 +  kModel(5)*dID.^5;

figure
p1 = axes;
plot(timeID(1:Fskip:end-1),consForce1(1:Fskip:end),'k')
hold on
plot(timeID(1:Fskip:end),ModFc(1:Fskip:end),'c--')
% plot(timeID(1:Fskip:end),model2(DispID(1:Fskip:end)),'g--')
tin = 0.5; tfin = 2.5; yin = -11; yfin = 11;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
xlabel('Time [s]')
ylabel('Force [N]')
xlim([0 45]);ylim([-12 12]);
% title("Experimental Conservative Force")
legend('Experiment','Identified Model', 'Interpreter', 'latex','location','southeast')

p2 = axes;
p2.Position = [0.40 0.63 0.45 0.25];
plot(timeID(1:Fskip:end-1),consForce1(1:Fskip:end),'k')
hold on
plot(timeID(1:Fskip:end),ModFc(1:Fskip:end),'c--')

xlim([tin tfin])
% ylim([-7 7])
annotation('Arrow','Position',[0.22679 0.75 0.11429 0], ...
    'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')


xMod = -0.01:1e-4:0.01;
model1 = K(1)*xMod + K(2)*xMod.^2 + K(3)*xMod.^3 + K(4)*xMod.^4 + K(5)*xMod.^5;                  % by CL

figure
% plot(DispID(Fstart:Fskip:end-1),consForce1(Fstart:Fskip:end),'k');
plot(DispID(1:Fskip:end-1),consForce1(1:Fskip:end),'k');
hold on
% plot(H(1:end-1,1),H(1:end-1,2),'r-')
% plot(xMod,model1(xMod),'r')
plot(xMod,model1,'r')
xlim([-0.008 0.008]);ylim([-13 13])
% plot(xMod,model2(xMod),'b-.')
xlabel('Displacement [m]')
ylabel('Force [N]')
legend('Experiment','Identified Model', 'Interpreter', 'latex','Location','northwest')
set(gcf,'renderer','painters')

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%%

IC = [0 0];
TF = time(1:sum(time <= 0.01));
Fsim1 = Force(1:sum(time <= 0.01));
Fsim2 = Data.Force(1:sum(time <= 0.01),3);
Fsim3 = Data.Force(1:sum(time <= 0.01),7);
Fsim4 = Data.Force(1:sum(time <= 0.01),20);
Fsim5 = Data.Force(1:sum(time <= 0.01),25);


[~,yCase1] = ode45(@(t,y) sys(t,y,mass,kModel,dampModel1,TF,Fsim1),time,IC);

[~,yCase2] = ode45(@(t,y) sys(t,y,mass,kModel,dampModel1,TF,Fsim2),time,IC);
[~,yCase3] = ode45(@(t,y) sys(t,y,mass,kModel,dampModel1,TF,Fsim3),time,IC);
[~,yCase4] = ode45(@(t,y) sys(t,y,mass,kModel,dampModel1,TF,Fsim4),time,IC);
[~,yCase5] = ode45(@(t,y) sys(t,y,mass,kModel,dampModel1,TF,Fsim5),time,IC);

%------------- Comparison

a_err = norm(yCase1(:,1)*1000-Disp*1000)/norm(Disp*1000);

MoDAL.PlotTSWTFT_Compare(time,Disp*1000,time,yCase1(:,1)*1000,0,20,"label", ...
    'Dispmm','title',['Displacement Response, Force = ' num2str(round(max(Fsim1))) ' N'], ...
    'fontSize',14,'legends',{'Experiment','Model - EDDI'},'timeEnd',45)
set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

MoDAL.PlotTSWTFT_Compare(time,Data.Disp(:,7)*1000,time,yCase3(:,1)*1000,0,20,"label", ...
    'Dispmm','title',['Displacement Response, Force = ' num2str(round(max(Fsim3))) ' N'], ...
    'fontSize',14,'legends',{'Experiment','Model'},'timeEnd',45)
set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

MoDAL.PlotTSWTFT_Compare(time,Data.Disp(:,25)*1000,time,yCase5(:,1)*1000,0,20,"label", ...
    'Dispmm','title',['Displacement Response, Force = ' num2str(round(max(Fsim5))) ' N'], ...
    'fontSize',14,'legends',{'Experiment','Model'},'timeEnd',45)

toc
set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

function [dy] = sys(t,y,m,k,c,TF,Force)
if t <= TF(end)
    F = interp1(TF,Force,t);
else 
    F = 0;
end

dy(1,1) = y(2);
dy(2,1) = -1/m*(c(1)*y(2) + c(2)*y(1)^2*y(2)+ c(3)*y(2)^2*y(2) + ...
                k(1)*y(1)+k(2)*y(1)^2+k(3)*y(1)^3+k(4)*y(1)^4+k(5)*y(1)^5-F);

end
