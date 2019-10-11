function [x,v,a]=NewmarkBeta_L_Step(M_all,K_all,C_all,F,x0,v0,dt,RecordLength)
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
% a(:,1) = inv(M_all(:,:,1))*(F(:,1)-C_all(:,:,1)*v0-K_all(:,:,1)*x0);
a(:,1) = inv(M_all(1))*(F(:,1)-C_all(1)*v0-K_all(1)*x0);

deta = 0.5;
alpha = 0.25;

a0 = 1/alpha/dt^2;
a1 = deta/alpha/dt;
a2 = 1/alpha/dt;
a3 = 1/2/alpha-1;
a4 = deta/alpha-1;
a5 = dt*(deta/alpha-2)/2;
a6 = 1-deta/alpha;
a7 = (1-deta/alpha/2)*dt;



for i=1:RecordLength-1
    %-------S-DOF---------
    M = M_all(i+1);
    K = K_all(i+1);
    C = C_all(i+1);
%-------M-DOF---------
%     M = M_all(:,:,i);
%     K = K_all(:,:,i);
%     C = C_all(:,:,i);
    K_f = K+a0*M+a1*C;
    iK = inv(K_f);
    
    F_f(:,i+1) = F(:,i+1)+M*(a0*x(:,i)+a2*v(:,i)+a3*(a(:,i)))+C*(a1*x(:,i)+a4*v(:,i)+a5*(a(:,i)));
    x(:,i+1) = iK*F_f(:,i+1);
    a(:,i+1) = a0*(x(:,i+1)-x(:,i))-a2*v(:,i)-a3*a(:,i);
    v(:,i+1) = a1*(x(:,i+1)-x(:,i))+a6*v(:,i)+a7*a(:,i);  
end
end