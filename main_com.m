clc;clear all;close all

% B=[1.2];
% f=[0.3];
% theta=[-pi/4];
% xi=[-0.001];
% %----------------
% dt=0.1;
% npo=200;
% t=[0:1:npo-1]*dt;
% %
% for nn=1:length(B)
%    Comp(nn,:)=B(nn)*exp(xi(nn)*t).*cos(2*pi*f(nn)*t+theta(nn));
% end
% sig0=sum(Comp,1);
% m=ones(length(t))*1;   % mass
% c=ones(length(t))*2;    % damping
% k=ones(length(t))*4;   % stiffness

% m=cos(2*t)+4;
% k=sin(5*t)+4;
% c=2*t+3;
% m=1;c=2;k=2;
% x0=0.2;
% v0=0.4;
m=[2,0;0,1];
k=[6,-2;-2,4];
c=[0,0;0,0];
dt=0.1;
npo=100;
t=[0:1:npo-1]*dt;
F=[0+t*0;10+t*0];
x0=[0;0];
v0=[0;1];
theta=1.4;
% FT=[sin(2*(t+theta*dt));3*cos(t+theta*dt)];
FE=@Fun;



[xhb,vhb,ahb]=Houbolt(m,k,c,F,x0,v0,dt,length(t));

% [xnl,vnl,anl]=NewmarkBeta_L(m,k,c,F,x0,v0,dt,length(t));
[xrk,vrk,ark]=RungeKuttaIV(m,k,c,FE,x0,v0,dt,length(t));

[xwt,vwt,awt]=Wilson_theta(m,k,c,FE,x0,v0,dt,length(t),theta);

% 
[xnl,vnl,anl]=NewmarkBeta_L(m,k,c,F,x0,v0,dt,length(t));
[xcd,vcd,acd]=CentralDifferenceM(m,k,c,F,dt,x0,v0,length(t));
[xnn,vnn,ann]=NewmarkBeta_NonL(m,k,c,F,dt,x0,v0,length(t));
%
% % [xnl,vnl,anl]=NewmarkBeta_L_Step(m,k,c,sig0,x0,v0,dt,length(t));
% % [xcd,vcd,acd]=CentralDifferenceM_Step(m,k,c,sig0,dt,x0,v0,length(t));
% % [xnn,vnn,ann]=NewmarkBeta_NonL_Step(m,k,c,sig0,dt,x0,v0,length(t));
% 
% figure(5)
% plot(t,xnl(1,:),'-g','linewidth',3);
% xlabel('Time(s)');ylabel('Displacement(m)');title('Newmark-beta Method')
% 
% 
% figure(7)
% plot(t,anl(1,:),'-m','linewidth',3);
% xlabel('Time(s)');ylabel('Accleration(m/s^2)');title('Newmark-beta Method')
% 
% 
% figure(6)
% plot(t,xcd(1,:),'-r','linewidth',3);
% xlabel('Time(s)');ylabel('Displacement(m)');title('Central Difference Method')
% 
% 
% figure(8)
% plot(t,acd(1,:),'-b','linewidth',3);
% xlabel('Time(s)');ylabel('Accleration(m/s^2)');title('Central Difference Method')
% 
xnl_1 = anl(1,:);
xcd_1 = acd(1,:);
xnn_1 = ann(1,:);
xwt_1 = awt(1,:);
xhb_1 = ahb(1,:);
xrk_1 = ark(1,:);

xnl_2 = anl(2,:);
xcd_2 = acd(2,:);
xnn_2 = ann(2,:);
xwt_2 = awt(2,:);
xhb_2 = ahb(2,:);
xrk_2 = ark(2,:);

figure(1)
plot(t,xnl_1,'-g','linewidth',3);
hold on
plot(t,xcd_1,'--r','linewidth',2.5);
hold on
plot(t,xnn_1,'-.m','linewidth',3.5);
hold on
plot(t,xwt_1,'-+b','linewidth',1.0);
hold on
plot(t,xhb_1,'-oc','linewidth',1.5);
hold on
plot(t,xrk_1,':sm','linewidth',1.5);


plot(t,xnl_2,'-g','linewidth',3);
hold on
plot(t,xcd_2,'--r','linewidth',2.5);
hold on
plot(t,xnn_2,'-.m','linewidth',3.5);
hold on
plot(t,xwt_2,'-+b','linewidth',1.0);
hold on
plot(t,xhb_2,'-oc','linewidth',1.5);
hold on
plot(t,xrk_2,':sm','linewidth',1.5);
l1=legend('Newmark Linear 1st','CentralDifference 1st','Newmark Non-Linear','Wilson-theta 1st','Houbolt 1st','Runge-Kutta 1st','Newmark Linear 2nd','CentralDifference 2nd','Newmark Non-LinearCJF 2nd','Wilson-theta 2nd','Houbolt 2nd','Runge-Kutta 2nd');
set(l1,'Fontname', 'Times New Roman','FontWeight','bold','FontSize',12)
% vnl_1 = vnl(1,:);
% vcd_1 = vcd(1,:);
% vnn_1 = vnn(1,:);
% 
% 
% figure
% plot(t,vnl_1,'-g','linewidth',3);
% hold on
% plot(t,vcd_1,'--r','linewidth',2.5);
% hold on
% plot(t,vnn_1,'-.m','linewidth',3.5);
% legend('Newmark Linear','CentralDifference','Newmark Non-Linear');

% anl_1 = anl(1,:);
% acd_1 = acd(1,:);
% ann_1 = ann(1,:);
% 
% figure
% plot(t,anl_1,'-g','linewidth',3);
% hold on
% plot(t,acd_1,'--r','linewidth',2.5);
% hold on
% plot(t,ann_1,'-.m','linewidth',3.5);
% legend('Newmark Linear','CentralDifference','Newmark Non-Linear');
% return