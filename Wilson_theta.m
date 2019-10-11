function [x,v,a]=Wilson_theta(M,K,C,FE,x0,v0,dt,RecordLength,theta)
%================Newmark-Beta (Linear system)method=====================
%           obtain the response of the dynamic system
%           [x,v,a]=newmarkb(M,K,C,N,F,x0,v0,a0,dt,RecordLength,abar)
%           M - mass matrix
%           K - stiffness matrix
%           C - damping matrix
%           F - loads
%           x0 - initial displacement
%           v0 - initial velocity
%           dt - interval
%           RecordLength - number of sampling points
%           theta - the value of the parameter of Wilson theta
t=(0:dt:(RecordLength-1)*dt);
F=Fun(t);
FT=Fun(t+theta*dt);

a0=6/((theta*dt)^2);
a1=3/(theta*dt);
a2=2*a1;
a3=theta*dt/2;
a4=a0/theta;
a5=-a2/theta;
a6=1-(3/theta);
a7=dt/2;
a8=dt^2/6;

Ks=K+a0*M+a1*C;

U=zeros(size(F,1),RecordLength);
V=zeros(size(F,1),RecordLength);
A=zeros(size(F,1),RecordLength);
FTs=zeros(size(F,1),RecordLength);
UT=zeros(size(F,1),RecordLength);
U(:,1)=x0;
V(:,1)=v0;
A(:,1) = inv(M)*(F(:,1)-C*v0-K*x0);
for i=1:RecordLength-1

    FTs(:,i)=F(:,i)+theta*(FT(:,i)-F(:,i))+M*(a0*U(:,i)+a2*V(:,i)+2*A(:,i))+C*(a1*U(:,i)+2*V(:,i)+a3*A(:,i));
    UT(:,i)=Ks\FTs(:,i);
    
    A(:,i+1)=a4*(UT(:,i)-U(:,i))+a5*V(:,i)+a6*A(:,i);
    V(:,i+1)=V(:,i)+a7*(A(:,i+1)+A(:,i));
    U(:,i+1)=U(:,i)+dt*V(:,i)+a8*(A(:,i+1)+2*A(:,i));
end

x=U;
v=V;
a=A;
end