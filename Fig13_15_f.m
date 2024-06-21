%% 
clc,clear,close all

load("SH_data.mat")

%%
clc,close all
[(1:length(Data.MaxForces))' Data.MaxForces]

Expi = 1; % 

time = Data.Time;
Acc = detrend(Data.Acc(:,Expi));

Fs = 1/(time(2)-time(1));

Vel = Data.Vel(:,Expi);
Disp = Data.Disp(:,Expi);
Force = Data.Force(:,Expi);
MaxF = Data.MaxForces(Expi);

%%

MoDAL.PlotForce(time,Force,"fontSize",14,"timeStart",0);
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

[~,pend] = min(abs(time-2));                                                                % Added by CL

% Find Time of Peak Velocity
[~,vMaxInd] = max(abs(Vel));

timeID = time(vMaxInd:pend); AccID = Acc(vMaxInd:pend); VelID = Vel(vMaxInd:pend);DispID = Disp(vMaxInd:pend);
time = time(vMaxInd:pend)-time(1);
hop = 1;
timeID = timeID(1:hop:end); AccID = AccID(1:hop:end); VelID = VelID(1:hop:end);DispID = DispID(1:hop:end);

MoDAL.PlotTSWTFT(timeID,1000*DispID,0,50,"label",'Dispmm','title', ...
    ['Displacement Response, Force = ' num2str(round(MaxF)) ' N'], ...
    "fontSize",14,"timeEnd",2,"numFreq",200)

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

[~,idx] = findpeaks(-abs(DispID));

mass = 0.088; 
KEId = 1/2*mass*VelID.^2;
R = KEId(1)-KEId(idx);

% Plot Kinetic Energy and KE when displacement is zero.
figure
plot(timeID,KEId,'k',timeID(idx),KEId(idx),'go')
% plot(timeID,KEId,'k',t0,KEVelD0,'go','MarkerSize',5);
hold on
tin = 0.1; tfin = 0.5; yin = 0.000; yfin = 0.004;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
xlim([0 2])
% title("Experimental Kinetic Energy")
xlabel('Time [s]'); ylabel('Energy [J]');
legend({'Experimental KE, T(t)', 'T(t=\gamma_i)'},Location='southeast'); % , 'Interpreter', 'latex',Location='southeast'
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

p2 = axes;
p2.Position = [0.40 0.35 0.45 0.5];
plot(timeID,KEId,'k'); hold on; plot(timeID(idx),KEId(idx),'go','MarkerSize',5);
xlim([tin tfin]); ylim([yin yfin])
annotation('Arrow','Position',[0.235,0.2263,0.0986,0.0920],'color',[0.7 0.7 0.7])
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

dampModel1 = Am1\R

% Compute Model Forces & Energies
dissForce1 = dampModel1(1)*VelID+dampModel1(2)*DispID.^2.*VelID+dampModel1(3)*VelID.^3;

dissEnergy1 = cumtrapz(timeID,VelID.*dissForce1);

MEEst_Model1 = KEId(1)-dissEnergy1;

% Plot dissipated energies
figure
plot(timeID(idx),R,'ko',timeID,dissEnergy1,'r')
hold on
tin = 0.1; tfin = 0.5; yin = 0.022; yfin = 0.028;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
xlim([0 2]); %ylim([0 0.014])
xlabel('Time [s]')
ylabel('Energy [J]')
% title("Comparison of Dissipated Energies")
legend('Experimental Dissipated Energy','Estimated Dissipated Energy','location', ...
    'southeast')
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')


p2 = axes;
p2.Position = [0.42 0.39048 0.45 0.42];
plot(timeID(idx),R,'ko',timeID,dissEnergy1,'r')
xlim([0.1 0.5])
h = annotation('Arrow','Position',[0.22679 0.708 0.10186 -0.12],'color',[0.7 0.7 0.7]); 
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')


%% Plot mechanical energies
figure
plot(timeID,KEId,'k',timeID,MEEst_Model1,'r')
hold on
tin = 0.1; tfin = 0.5; yin = 0; yfin = 0.0045;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
xlim([0 2])
ylim([0 0.03])
xlabel('Time [s]');ylabel('Energy [J]')
% title("Comparison of Identified Mechanical Energies")
% title(['Force = ' num2str(MaxF) ' N'])
legend('Experimental KE','Estimated ME', 'Interpreter', 'latex','position', ...
    [0.50009 0.13801 0.34732 0.15595])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')


p2 = axes;
p2.Position = [0.40 0.36 0.45 0.5];
plot(timeID,KEId,'k',timeID,MEEst_Model1,'r')
xlim([0.1 0.5]);
annotation('Arrow','Position',[0.22679 0.24143 0.12143 0.16667], ...
    'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')


%%
Fstart = 1;
Fskip = 1;
potEnergy1 = MEEst_Model1 - KEId;

% Plot potential energy
figure
plot(timeID,potEnergy1,'k')
hold on
tin = 0.1; tfin = 0.5; yin = -0.002; yfin = 0.006;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
legend('Estimated Potential Energy','position',[0.41429 0.1369 0.44643 0.045238]);
xlim([0 4]);ylim([-.003 .014])
xlabel('Time [s]'); ylabel('Energy [J]'); 
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')


p2 = axes;
p2.Position = [0.425 0.375 0.45 0.5];
plot(timeID,potEnergy1,'k')
xlim([tin tfin]);ylim([yin .006])
h = annotation('Arrow','Position',[0.235 0.39762 0.13464 0.028572], ...
    'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%%
smoothValue = 100;

L = 2*KEId-MEEst_Model1;
dL_dx = gradient(L)./gradient(DispID);


m_a =  mass*AccID;  

consForce1 = -(dL_dx-m_a);

consForce1 = detrend(smoothdata(consForce1,'movmedian',smoothValue));     % Modified by CL
%% Identify model for conservative force

dID = DispID(Fstart:end);

x = dID;
theta = [x x.^2 x.^3 x.^4 x.^5];

K = theta\(consForce1)

kModel1 = K;

% kModel1 = [model1.k1; model1.k2; model1.k3; model1.k4; model1.k5]
ModFc = kModel1(1)*dID'+kModel1(2)*(dID').^2+kModel1(3)*(dID').^3+kModel1(4)*(dID').^4+kModel1(5)*(dID').^5;

figure
p1 = axes;
plot(timeID(Fstart:Fskip:end),consForce1(Fstart:Fskip:end),'k')
hold on
plot(timeID(Fstart:Fskip:end),ModFc(1:Fskip:end),'c--')
% plot(timeID(1:Fskip:end),model2(DispID(1:Fskip:end)),'g--')
tin = 0.1; tfin = 0.5; yin = -4; yfin = 4;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
xlabel('Time [s]')
ylabel('Force [N]')
xlim([0 2]);ylim([-6 9]);
% title("Experimental Conservative Force")
legend('Experiment','Identified Model', 'Interpreter', 'latex','location','southeast')
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

p2 = axes;
p2.Position = [0.40 0.6 0.45 0.25];
plot(timeID(1:Fskip:end),consForce1,'k')
hold on
plot(timeID(Fstart:Fskip:end),ModFc,'c--')
set(gca,'fontsize',12)
xlim([tin tfin]); ylim([yin yfin])
annotation('Arrow','Position',[0.22679 0.656 0.11429 0.07], ...
    'color',[0.7 0.7 0.7])

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

xMod = -0.01:1e-4:0.01;
model1 = K(1)*xMod + K(2)*xMod.^2 + K(3)*xMod.^3+ K(4)*xMod.^4 + K(5)*xMod.^5;

%%
xi = 1;
figure
% plot(DispID(Fstart:Fskip:end-1),consForce1(Fstart:Fskip:end),'k');
plot(DispID(xi:Fskip:end),consForce1(xi:Fskip:end),'k');
hold on
% plot(H(1:end-1,1),H(1:end-1,2),'r-')
plot(xMod,model1,'r','LineWidth',2)
xlim([-0.008 0.008]);ylim([-10 10])
% plot(xMod,model2(xMod),'b-.')
xlabel('Displacement [m]')
ylabel('Force [N]')
legend('Experiment','Identified Model', 'Interpreter', 'latex','Location','northwest')
set(gca,'fontsize',12)
set(gcf,'renderer','painters')
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%%

IC = [0 0];
TF = time(1:sum(time <= 0.01));
% close all
S2 = 2;
S4 = 3;

Fsim1 = Force(1:sum(time <= 0.01));
Fsim2 = Data.Force(1:sum(time <= 0.01),S2);
Fsim4 = Data.Force(1:sum(time <= 0.01),S4);

Fsim1 = Fsim1-Fsim1(1);
Fsim2 = Fsim2-Fsim2(1);
Fsim4 = Fsim4-Fsim4(1);

[~,yCase1] = ode45(@(t,y) sys2(t,y,mass,kModel1,dampModel1,TF,Fsim1),time,IC);
[~,yCase2] = ode45(@(t,y) sys2(t,y,mass,kModel1,dampModel1,TF,Fsim2),time,IC);
[~,yCase4] = ode45(@(t,y) sys2(t,y,mass,kModel1,dampModel1,TF,Fsim4),time,IC);

Disp = Disp(1:length(time));
Disp2 = Data.Disp(:,S2);     Disp2 = Disp2(1:length(time));
Disp3 = Data.Disp(:,S4);     Disp3 = Disp3(1:length(time));

MoDAL.PlotTSWTFT_Compare(time,Disp*1000,time,yCase1(:,1)*1000,0,50,"label", ...
    'Dispmm','title',['Displacement Response, Force = ' num2str(round(max(Fsim1))) ' N'], ...
    'fontSize',14,'legends',{'Experiment','Model - EDDI'},'timeEnd',2)
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

MoDAL.PlotTSWTFT_Compare(time,Disp2*1000,time,yCase2(:,1)*1000,0,50,"label", ...
    'Dispmm','title',['Displacement Response, Force = ' num2str(round(max(Fsim2))) ' N'], ...
    'fontSize',14,'legends',{'Experiment','Model'},'timeEnd',2)
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

MoDAL.PlotTSWTFT_Compare(time,Disp3*1000,time,yCase4(:,1)*1000,0,50,"label", ...
    'Dispmm','title',['Displacement Response, Force = ' num2str(round(max(Fsim4))) ' N'], ...
    'fontSize',14,'legends',{'Experiment','Model'},'timeEnd',2)
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%%

function dy = sys(t,y,m,k,c,TF,Force)
if t <= TF(end)
    F = interp1(TF,Force,t);
else 
    F = 0;
end

dy(1,1) = y(2);
dy(2,1) = -1/m*(c(1)*y(2) +  c(2)*y(1)^2*y(2)+ c(3)*y(2)^2*y(2) + ...
                k(1)*y(1)+k(2)*y(1)^3+k(3)*y(1)^5-F);

end


function dy = sys2(t,y,m,k,c,TF,Force)
if t <= TF(end)
    F = interp1(TF,Force,t);
else 
    F = 0;
end

dy(1,1) = y(2);
dy(2,1) = -1/m*(c(1)*y(2) +  c(2)*y(1)^2*y(2)+ c(3)*y(2)^2*y(2) + ...
                k(1)*y(1)+k(2)*y(1)^2+k(3)*y(1)^3+k(4)*y(1)^4+k(5)*y(1)^5-F);

end
