function [x v a]=NewmarkBeta_NonL(M,K,C,F,dt,x0,v0,RecordLength)
%==============Newmark-Beta (Non-Linear system) method=====================
%           obtain the response of the dynamic system
%           [x,v,a]=NewmarkBeta_NonL(M,K,C,F,dt,x0,v0,a0,RecordLength,abar)
%           M - mass matrix
%           K - stiffness matrix
%           C - damping matrix
%           F - loads
%           x0 - initial displacement
%           v0 - initial velocity
%           dt - interval
%           RecordLength - number of sampling points
%           abar - 
%==========================================================================
delta_t = dt;

x(:,1) = x0;
v(:,1) = v0;
% a0 = inv(M)*(F-C*v0-K*x0);
a(:,1) = inv(M)*(F(:,1)-C*v0-K*x0);
gamma = 0.5;
beta = 0.25;

am1 = 1/beta/delta_t^2;
am2 = 1/beta/delta_t;
am3 = 1/2/beta;

ac1 = gamma/beta/delta_t;
ac2 = gamma/beta;
ac3 = (gamma/2/beta-1)*delta_t;


K_tr = K+am1*M+ac1*C;
iK = inv(K_tr);
for i=1:RecordLength-1
    F_delta(:,i) = F(:,i+1)-F(:,i);
%     F_f(:,i) = F_delta(:,i)+M*(am2*v(:,i)+am3*(a(:,i)+abar))+C*(ac2*v(:,i)+ac3*(a(:,i)+abar));
    F_f(:,i) = F_delta(:,i)+M*(am2*v(:,i)+am3*(a(:,i)))+C*(ac2*v(:,i)+ac3*(a(:,i)));
    x_delta(:,i) = iK*F_f(:,i);
    v_delta(:,i) = ac1*x_delta(:,i)-ac2*v(:,i)-ac3*a(:,i);
    a_delta(:,i) = am1*x_delta(:,i)-am2*v(:,i)-am3*a(:,i);
    x(:,i+1) = x(:,i)+x_delta(:,i);
    v(:,i+1) = v(:,i)+v_delta(:,i);
    a(:,i+1) = a(:,i)+a_delta(:,i);
end
end