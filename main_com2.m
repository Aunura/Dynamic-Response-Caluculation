clc;clear all;close all

B=[1.2];
f=[0.3];
theta=[-pi/4];
xi=[-0.01];
%----------------
dt=0.01;
npo=2000;
t=[0:1:npo-1]*dt;
%
for nn=1:length(B)
   Comp(nn,:)=B(nn)*exp(xi(nn)*t).*cos(2*pi*f(nn)*t+theta(nn));
end
sig0=sum(Comp,1);
%----------------
x0=0;v0=0;
rand('state',0);
val_variant_m=0.1;
val_variant_k=0.35;
val_variant_c=0.25;
m=4;k=2;c=3;
m_Time=m+rand(1,npo).*val_variant_m
k_Time=k+rand(1,npo).*val_variant_k
c_Time=c+rand(1,npo).*val_variant_c

% m=cos(2*t)+4;
% k=sin(5*t)+4;
% c=2*t+3;

x0=0;v0=0;

% [xnl,vnl,anl]=NewmarkBeta_L(m,k,c,sig0,x0,v0,dt,length(t));
% [xcd,vcd,acd]=CentralDifferenceM(m,k,c,sig0,dt,x0,v0,length(t));
% [xnn,vnn,ann]=NewmarkBeta_NonL(m,k,c,sig0,dt,x0,v0,length(t));

[xnl,vnl,anl]=NewmarkBeta_L_Step(m_Time,k_Time,c_Time,sig0,x0,v0,dt,length(t));
[xcd,vcd,acd]=CentralDifferenceM_Step(m_Time,k_Time,c_Time,sig0,dt,x0,v0,length(t));
[xnn,vnn,ann]=NewmarkBeta_NonL_Step(m_Time,k_Time,c_Time,sig0,dt,x0,v0,length(t));


xnl_1 = xnl(1,:);
xcd_1 = xcd(1,:);
xnn_1 = xnn(1,:);

figure
plot(t,xnl_1,'-g','linewidth',3);
hold on
plot(t,xcd_1,'--r','linewidth',2.5);
hold on
plot(t,xnn_1,'-.m','linewidth',3.5);
legend('Newmark Linear','CentralDifference','Newmark Non-Linear');


return