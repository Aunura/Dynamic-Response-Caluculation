function [x,v,a]=NewmarkBeta_L(M,K,C,F,x0,v0,dt,RecordLength)
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
%           abar - 
%=======================================================================
x(:,1) = x0;
v(:,1) = v0;
a(:,1) = inv(M)*(F(:,1)-C*v0-K*x0);

deta = 0.5;
alpha = 0.25;

for i=1:RecordLength-1
%     F_f(:,i+1)=F(:,i+1)+M*(a0*x(:,i)+a2*v(:,i)+a3*(a(:,i)+abar))+C*(a1*x(:,i)+a4*v(:,i)+a5*(a(:,i)+abar));
    F_f(:,i+1) = F(:,i+1)+M*(a0*x(:,i)+a2*v(:,i)+a3*(a(:,i)))+C*(a1*x(:,i)+a4*v(:,i)+a5*(a(:,i)));
    x(:,i+1) = iK*F_f(:,i+1);
    a(:,i+1) = a0*(x(:,i+1)-x(:,i))-a2*v(:,i)-a3*a(:,i);
    v(:,i+1) = a1*(x(:,i+1)-x(:,i))+a6*v(:,i)+a7*a(:,i);  
end
end