clc,clear,close all

% Parameters
mass = 0.05; d = 0.5; dnl = 4e3; k = 300; knl = 3e8; 

% Initial conditions
IC = [0 10];

dt = 1e-4;
t = (0:dt:1)';

Options = odeset; Options.RelTol = 1e-12; Options.AbsTol = 1e-16;
[~,y] = ode45(@(t,y) sys2(t,y,mass,k,d,dnl,knl),t,IC,Options);

MoDAL.PlotTSWTFT(t,y(:,1),0,200,"label",'Disp', ...
    "fontSize",14,"timeEnd",1.0)
set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')
%%
% Find times where displacement/velocity are minimized
[~,idx] = findpeaks(-abs(y(:,1)));
t = t(idx(1):end);
y = y(idx(1):end,:);
% y(1,:)
[To,idx] = findpeaks(-abs(y(:,1)));

T = 1/2*mass*y(:,2).^2;
V = 1/2*k*y(:,1).^2 + 1/4*knl*y(:,1).^4;
EM = T + V;

% ac = gradient(y(:,2))/dt;
% dEM = mass*y(:,2).*ac + k*y(:,1).*y(:,2) + knl*y(:,1).^3.*y(:,2);

% figure;
% plot(EM);hold on;plot(T,'k');plot(dEM,'b');plot(-(1./(dnl*y(:,2).^2 + d)).*dEM,'--r')
% legend('Mechanical Energy (ME)','Kinetic Energy','Derivative of the ME','Theoretical Derivative of the ME')

%............... LHS ...............
Pow = 3;
A = zeros(length(idx),Pow);
for m = 2 : Pow + 1
    xdot = y(:,2).^ m;
    At = cumtrapz(t,xdot);
    A(:,m-1) = At(idx);
end
x2_xdot = y(:,1).^ 2.*y(:,2).^2;
RR1 = cumtrapz(t,x2_xdot);
A1 = RR1(idx);

% xdot2_xdot = y(:,2).^ 2.*y(:,2).^2;
% RR2 = cumtrapz(t,xdot2_xdot);
% A2 = RR2(idx);

Am = [A A1];

%............... RHS ...............
R = T(1)-T(idx);

%..............   Damping model calculation  ...............
c = Am\R

%..............   Dissipative force identification  ...............

x = y(:,1);
v = y(:,2);

dissForce = c(1)*v + c(2)*v.^2 + c(3)*v.^3 + ...
             c(4)*x.^2.*v;

%% Plot linear damping and nonlinear damping forces from the simulation
figure
plot(t,c(1)*v,'k',t,c(4)*x.^2.*v,'r')
hold on
tin = .05; tfin = .15; yin = -3; yfin = 3;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
xlabel('Time [s]')
ylabel('Force [N]')
% title("Comparison of Dissipated Energies")
legend('Linear damping Force','Nonlinear damping Force','location','southeast')

p2 = axes;
p2.Position = [0.495 0.615 0.37 0.30];
plot(t,c(1)*v,'k',t,c(4)*x.^2.*v,'r')
xlim([tin tfin]);ylim([yin yfin])
annotation('Arrow','Position',[0.245 0.72 0.092857 0],'color',[0.7 0.7 0.7]) 
set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

%% Plot T, ME, T(gamma)
figure;
p1 = axes;
plot(t,T,'k'); hold on; plot(t,EM,'b');plot(t(idx),T(idx),'go','MarkerSize',5);
ylim([0 2.5])
xlabel('Time /s'); ylabel('Energy [J]');
tin = 0.075; tfin = 0.1; yin = 0.45; yfin = 0.7;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
p2 = axes;
p2.Position = [0.40 0.35 0.45 0.5];
plot(t,T,'k'); hold on; plot(t,EM,'b');plot(t(idx),T(idx),'go','MarkerSize',5);
xlim([tin tfin]);ylim([yin yfin])
legend('T', 'ME', 'T (\gamma_i)'); % set (legend_handle,'Interpreter','latex')
annotation('Arrow','Position',[0.21,0.345,0.08,0.10],'color',[0.7 0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

%% Plot Dissipative Energy
%..............   Dissipated energy  ...............
dissEnergy = cumtrapz(t,v.*dissForce);

figure
Real_DE = cumtrapz(t,v.*(d*y(:,2) + dnl.*y(:,1).^2.*y(:,2)));
plot(t,Real_DE ,'k',t,dissEnergy,'r--')
hold on
yin = 1.6; yfin = 1.9;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
ylim([0 2.5])
xlabel('Time [s]')
ylabel('Dissipative Energy [J]')
% title("Comparison of Dissipated Energies")
legend('Exact','Estimated','location', ...
    'southeast')

p2 = axes;
p2.Position = [0.46 0.38 0.38 0.42];
plot(t,Real_DE ,'k',t,dissEnergy,'r--')
xlim([tin tfin])
annotation('Arrow','Position',[.206 0.65 0.092857 -0.02619],'color',[0.7 0.7 0.7]) 
set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

%%  Plot Exact and Estimated ME
% ..................   ME calculation ......................
MEEst = T(1)-dissEnergy;
figure;
p1 = axes;
plot(t,EM,'b');hold on; plot(t,MEEst,'r--');
ylim([0 2.5])
xlabel('Time [s]'); ylabel('Energy [J]');
set(gca,'fontsize',14)
% tin = 10; tfin = 16; yin = 0.4; yfin = 0.8;
yin = 0.45; yfin = 0.7;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
p2 = axes;
p2.Position = [0.40 0.35 0.45 0.5];
plot(t,EM,'b');hold on; plot(t,MEEst,'r--');
xlim([tin tfin]);ylim([yin yfin])
legend('Exact','Estimated'); % set (legend_handle,'Interpreter','latex')
annotation('Arrow','Position',[0.21,0.345,0.08,0.10],'color',[0.7 0.7 0.7])

set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

%% F1Est = detrend(smoothdata(F1Est,'movmedian',ll));     % Modified by CL

V1 = MEEst - T;
% figure
% plot(t,V,'k',t,V1,'c--'); xlabel('Time [s]'); ylabel('Potential Energy [J]')
% ylim([-0.1 4.5])
% % title("Comparison of Dissipated Energies")
% legend('Simulated','Calculated ','location','northeast')

%------ Lagrangian
L = 2*T-MEEst;
acc = gradient(v)/dt;
F1Est = gradient(L)./gradient(x) - mass*acc;      
F1Est = -F1Est;

%..............    Proposed model  ...............
dID = x;

% Using least-squares
clc
theta = [x  x.^2  x.^3  x.^4  x.^5];
model1 = theta\F1Est;
%..............    Sparse regression ..............   
kModel1 = [model1(1); model1(2); model1(3); model1(4); model1(5)]

Freal = k*dID + knl*dID.^3;

z = linspace(-0.013,0.013,1000);
Fmodel = kModel1(1)*z + kModel1(2)*z.^2 + kModel1(3)*z.^3 + kModel1(4)*z.^4 + kModel1(5)*z.^5; % ;
    
figure; plot(dID,Freal,'k')
xlabel('Displacement [m]');ylabel('Force [N]')
hold on;plot(z,Fmodel,'r--'); xlim([-.015 .015]); ylim([-750 750])
legend('Exact','Model','Location','Best')
set(findall(gcf,'-property','FontSize'),'FontSize',18, 'FontName', 'Times New Roman')

IC = [x(1) v(1)];
    
%..............  Model calculation ...............
[~,ymodel] = ode45(@(t,y) sys(t,y,mass,kModel1,c),t,IC,Options);

% Fd = c(1)*ymodel(:,2);
% figure
% plot(t,Fd,'r--');xlabel('Time [s]'); ylabel('Displacement [m]')
% legend('Damping force','Location','Best')

%..............  Comparison ...............
a_err = norm(ymodel(:,1)-y(:,1))/norm(y(:,1));

MoDAL.PlotTSWTFT_Compare(t,y(:,1),t,ymodel(:,1),0,200,"label",'Disp', ...
    'fontSize',14,'legends',{'Exact System','Identified Model'},'timeEnd',1)

set(findall(gcf,'-property','FontSize'),'FontSize',16, 'FontName', 'Times New Roman')

((0.5-0.49999)/.5)*100
((4000-4000.2)/4000)*100
((300-317.32)/300)*100
((3e8-3.003e8)/3e8)*100

function dy = sys2(t,y,m,k,d,dnl,knl)
dy(1,1) = y(2);
dy(2,1) = -1/m*(d*y(2) + dnl.*y(1).^2.*y(2) + k*y(1) + knl*y(1)^3);
end

function dy = sys(t,y,m,Xi,b)
dy(1,1) = y(2);
dy(2,1) = -1/m*(b(1)*y(2) + b(2)*y(2)*y(2) + b(3)*y(2)^2*y(2) + ...
                b(4)*y(1)^2*y(2) +...
               Xi(1)*y(1) + Xi(2)*y(1)^2 + Xi(3)*y(1)^3 + Xi(4)*y(1)^4 + Xi(5)*y(1)^5);

end