clear;close all;clc;

% Parameters
mass = 2;   % mass of the pendulum bob
l = .8;      % length of the pendulum
b = 0.1;    % damping coefficient
g = 9.81;   % acceleration due to gravity

% Initial conditions
theta0 = pi/2;  % initial angle
theta_dot0 = 0; % initial angular velocity

% Time span
dt = 1e-2;
t = (0:dt:100)';

% Solve ODE
Options = odeset; Options.RelTol = 1e-12; Options.AbsTol = 1e-16;
[t, y] = ode45(@(t, y) pendulum_ode(t, y, b, mass, g, l), t, [theta0 theta_dot0],Options);

MoDAL.PlotTSWTFT(t,y(:,1),0,1,"label",'Theta', ...
    "fontSize",14,"timeEnd",100)

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

% Plot results
% figure;
% plot(t, y(:, 1));
% title('Damped Pendulum Simulation');
% xlabel('Time (s)');
% ylabel('Angular Displacement (rad)');

% Find times where displacement/velocity are minimized
tog = t;
yog = y;
[~,idx] = findpeaks(-abs(y(:,1)));
t = t(idx(1):end);
y = y(idx(1):end,:);

[To,idx] = findpeaks(-abs(y(:,1)));

T = 1/2 * mass * (l *y(:,2)).^2;
V = mass * g * l * (1-cos(y(:, 1)));
EM = T + V;

% figure;
% plot(t,EM);hold on;plot(t,T,'k');
% legend('Mechanical Energy (ME)','Kinetic Energy')

%............... LHS ...............
Pow = 3;
A = zeros(length(idx),Pow);
for m = 2 : Pow + 1
    xdot = (y(:,2)).^ m;
    At = cumtrapz(t,xdot);
    A(:,m-1) = At(idx);
end
x2_xdot = y(:,1).^ 2.*(y(:,2)).^2;
RR1 = cumtrapz(t,x2_xdot);
A1 = RR1(idx);

Am = [A A1];

%............... RHS ...............
R = T(1)-T(idx);

%..............   Damping model calculation  ...............
c = Am\R

%..............   Restoring force identification  ...............
x = y(:,1);
v = y(:,2);                                                              % For polar coord

MEEst = T(1)-cumtrapz(t,c(1)*v.^2 + c(2)*v.^3 + c(3)*v.^4 + c(4)*x.^2.*v.^2);

figure;
p1 = axes;
plot(t,T,'k'); hold on; plot(t,EM,'b');plot(t(idx),T(idx),'go','MarkerSize',5);
% ylim([0 3.2])
xlabel('Time [s]'); ylabel('Energy [J]')
tin = 20.5; tfin = 25; yin = 4.5; yfin = 6;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
p2 = axes;
p2.Position = [0.42 0.4 0.45 0.45];
plot(t,T,'k'); hold on; plot(t,EM,'b');plot(t(idx),T(idx),'go','MarkerSize',5);
xlim([tin tfin]);ylim([yin yfin])
legend('T', 'ME', 'T(\gamma_i)'); % set (legend_handle,'Interpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

annotation('Arrow','Position',[0.28,0.45,0.1,0.04],'color',[0.7 0.7 0.7])

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%% Plot linear damping force from the simulation
figure
plot(t,b*v,'k')
xlabel('Time [s]');ylabel('Force [N-m]')
% title("Comparison of Dissipated Energies")
legend('Linear damping Force','location','southeast')
% ylim([-0.7 0.7])
set(gca,'fontsize',12)

%%
figure;
p1 = axes;
plot(t,EM,'b');hold on; plot(t,MEEst,'r--');
% ylim([0 3.2])
xlabel('Time [s]'); ylabel('Energy [J]')
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
p2 = axes;
p2.Position = [0.42 0.4 0.45 0.45];
plot(t,EM,'b');hold on; plot(t,MEEst,'r--');
xlim([tin tfin]);ylim([yin yfin])
legend('Exact', 'Estimated'); % set (legend_handle,'Interpreter','latex')
annotation('Arrow','Position',[0.28,0.45,0.1,0.04],'color',[0.7 0.7 0.7])

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%% Dissipative Energy

% Compute Model Forces & Energies
dissForce1 = c(1)*v + c(2)*v.^2 + c(3)*v.^3 + c(4)*x.^2.*v;

dissEnergy1 = cumtrapz(t,v.*dissForce1);

% MEEst_Model1 = T(1)-dissEnergy1;

% Plot dissipated energies
figure
Real_DE = cumtrapz(t,v.*(b*l^2*y(:,2)));
plot(t,Real_DE ,'k',t,dissEnergy1,'r--')
hold on
yin = 9.4; yfin = 10.6;
plot([tin tfin],[yin yin],'color',[0.7 0.7 0.7])
plot([tin tfin],[yfin yfin],'color',[0.7 0.7 0.7])
plot([tin tin],[yin yfin],'color',[0.7 0.7 0.7])
plot([tfin tfin],[yin yfin],'color',[0.7 0.7 0.7])
% ylim([0 3.2])
xlabel('Time [s]')
ylabel('Dissipative Energy [J]')
% title("Comparison of Dissipated Energies")
legend('Exact','Estimated','location', ...
    'southeast')
set(gca,'fontsize',12)

p2 = axes;
p2.Position = [0.5375 0.35048 0.38 0.42];
plot(t,Real_DE ,'k',t,dissEnergy1,'r--')
xlim([tin tfin])
annotation('Arrow','Position',[0.234 0.688 0.092857 -0.02619],'color',[0.7 0.7 0.7]) 

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%% Plot linear damping force from the simulation
figure
plot(t,b*v,'k')
xlabel('Time [s]');ylabel('Force [N-m]')
% title("Comparison of Dissipated Energies")
legend('Linear damping Force','location','southeast')
% ylim([-0.7 0.7])
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%% F1Est = detrend(smoothdata(F1Est,'movmedian',ll));     % Modified by CL

V1 = MEEst - T;

figure
plot(t,V,'k',t,V1,'c--'); xlabel('Time [s]'); ylabel('Energy [J]')
% ylim([-0.1 2.8])
% title("Comparison of Dissipated Energies")
legend('Real V','Calculated V','location','northeast')
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

L = 2*T-MEEst;
acc = gradient(v)/dt;
F1Est = gradient(L)./gradient(x) - mass*l^2*acc;                             
F1Est = -detrend(smoothdata(F1Est,'movmedian',20));         % Modified by CL

% figure;plot(F1Est)

%..............    Proposed model  ...............
% F1Est_real = (mass*g*l)*sin(x);F1Est_real = F1Est_real;figure;plot(F1Est_real);hold on;plot(F1Est,'r--')
% _real(1:end-1)
dID = x;

anonmodel1 = @(k1,k2,k3,x) k1*sin(x)+k2*cos(x)+k3*tan(x);
model1 = fit(dID,F1Est,anonmodel1, ...
    'StartPoint',[.5 .5 .5]);

anonmodel1 = @(k1,k2,k3,k4,k5,x) k1*x+k2*x.^2+k3*x.^3+k4*x.^4+k5*x.^5;
model2 = fit(dID,F1Est,anonmodel1, ...
    'StartPoint',[.5 .5 .5 .5 .5]);

theta = [x  x.^2  x.^3  x.^4  x.^5];
K = theta\F1Est

dID = x;
anonmodel1 = @(k1,k3,k5,x) k1*x+k3*x.^3+k5*x.^5;
model3 = fit(dID,F1Est,anonmodel1, ...
    'StartPoint',[1 1 1])


Stiffness = {'k1 = mgl','k2 = 0','k3 = -mgl/3!','k4 = 0','k5 = mgl/5!'}';
Exact = mass*g*l*[1;0;-1./factorial(3);0;1./factorial(5)];
kModel2 = [model2.k1; model2.k2; model2.k3; model2.k4; model2.k5];
% kModel3 = [model3.k1; 0; model3.k3; model3.k5]

T2 = table(Stiffness,Exact,K,kModel2)


%..............  Model Comparison  ..............
Stiffness = {'k1 = mgl, sin(theta)','k2 = 0, cos(theta)','k3 = 0, tan(theta)'}';
Exact = [mass*g*l;0;0];
kModel1 = [model1.k1; model1.k2; model1.k3];
T1 = table(Stiffness,Exact,kModel1)

Stiffness = {'k1 = mgl','k2 = 0','k3 = -mgl/3!','k4 = 0','k5 = mgl/5!'}';
Exact = mass*g*l*[1;0;-1./factorial(3);0;1./factorial(5)];
kModel2 = [model2.k1; model2.k2; model2.k3; model2.k4; model2.k5];
% kModel3 = [model3.k1; 0; model3.k3; model3.k5]

T2 = table(Stiffness,Exact,kModel2)


%..............  Force Comparison  ..............
z = linspace(-max(abs(y(:,1))),max(abs(y(:,1))),1000);
Fexact = mass*g*l*sin(z);
Fmod = kModel1(1)*sin(z) + kModel1(2)*cos(z) + kModel1(3)*tan(z); % ;
% 
% figure; plot(z,Fexact,'k',z,Fmod,'r')
% xlabel('$x$', 'interp','Latex');ylabel('Restoring force [N]')
% 
% legend('Exact','Model','Location','Best');

Fexact = mass*g*l*(z-z.^3/factorial(3)+z.^5/factorial(5));
Fmod = kModel2(1)*z+kModel2(2)*z.^2+kModel2(3)*z.^3+kModel2(4)*z.^4+kModel2(5)*z.^5;

figure; plot(z,Fexact,'k',z,Fmod,'r')
xlabel('$x$', 'interp','Latex');ylabel('Restoring force [N]')
legend('Exact','Model','Location','Best');
    
set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

IC = yog(1,:);
    
%..............  Model calculation ...............
[t,ymodel] = ode45(@(t,y) sys(t,y,mass,kModel2,c,l),tog,IC,Options);

% Fd = c(1)*ymodel(:,2);
% figure
% plot(t,Fd,'r--');xlabel('Time [s]'); ylabel('Displacement [m]')
% legend('Damping force','Location','Best')

%..............  Comparison ...............

MoDAL.PlotTSWTFT_Compare(t,yog(:,1),t,ymodel(:,1),0,1,"label",'DispND', ...
    'fontSize',14,'legends',{'Simulation','Model'},'timeEnd',100)

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

((0.064006-0.064)/0.064)*100
((15.696-15.652)/15.7)*100
((-2.616+2.58)/2.61)*100
((0.1308-0.112)/0.13)*100

function dydt = pendulum_ode(t, y, b, mass, g, l)
dydt = zeros(2, 1);
dydt(1) = y(2);
dydt(2) = -(b/mass)*y(2) - (g/l)*sin(y(1));
end


function dy = sys(t,y,m,Xi,b,l)
dy(1,1) = y(2);
dy(2,1) = -1/(m*l^2)*(b(1)*y(2) + b(2)*y(2)^2 + b(3)*y(2)^3 + ...
                b(4)*y(1)^2*y(2) + Xi(1)*y(1) + Xi(2)*y(1)^2 + Xi(3)*y(1)^3 + ...
          Xi(4)*y(1)^4 + Xi(5)*y(1)^5);

end