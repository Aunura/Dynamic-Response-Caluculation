function [x v a]=CentralDifferenceM_Step(M_all,K_all,C_all,F,dt,x0,v0,RecordLength)
%=====================CentralDifferenceM method=======================
%       obtain the response of the dynamic system
%       [x,v,a]=CentralDifferenceM(M,K,C,F,dt,x0,v0,a0,RecordLength,abar)
%       M - mass matrix
%       K - stiffness matrix
%       C - damping matrix
%       F - loads
%       x0 - initial displacement
%       v0 - initial velocity
%       dt - interval
%       RecordLength - number of sampling points
%       abar - 
%================================================================
delta_t = dt;
x(:,1) = x0;
v(:,1) = v0;
% a(:,1) = inv(M_all(:,:,1))*(F(:,1)-C_all(:,:,1)*v0-K_all(:,:,1)*x0);
a(:,1) = inv(M_all(1))*(F(:,1)-C_all(1)*v0-K_all(1)*x0);
a0 = a(:,1);

K_tr = M_all(1)/(delta_t^2)+C_all(1)/2/delta_t;
% K_tr = M_all(:,:,1)/(delta_t^2)+C_all(:,:,1)/2/delta_t;
iK=inv(K_tr);

%=============解决起步问题（前二个时间点的求解）===================
x_neg1 = x0-v0*delta_t+(a0/2)*delta_t^2;
% F_tr(:,1) = F(:,1)-(M_all(:,:,1)/(delta_t^2)-C_all(:,:,1)/2/delta_t)*x_neg1-(K_all(:,:,1)-2*M_all(:,:,1)/(delta_t^2))*x(:,1);
F_tr(:,1) = F(:,1)-(M_all(1)/(delta_t^2)-C_all(1)/2/delta_t)*x_neg1-(K_all(1)-2*M_all(1)/(delta_t^2))*x(:,1);
x(:,2) = iK*F_tr(:,1);
v(:,2) = (x(:,2)-x(:,1))/2/delta_t;
a(:,2) = (x(:,2)-2*x(:,1)+x_neg1)/(delta_t^2);
%===================end========================================

for i=2:RecordLength-1
%-------S-DOF---------
    M = M_all(i+1);
    K = K_all(i+1);
    C = C_all(i+1);
%-------M-DOF---------
%     M = M_all(:,:,i);
%     K = K_all(:,:,i);
%     C = C_all(:,:,i);
    K_tr = M/(delta_t^2)+C/2/delta_t;
    iK=inv(K_tr);
    
    F_tr(:,i) = F(:,i)-(M/(delta_t^2)-C/2/delta_t)*x(:,i-1)-(K-2*M/(delta_t^2))*x(:,i);
%     F_tr(:,i) = F(:,i)-K*(x(:,i)+abar)+M*(2*x(:,i)/(delta_t^2)-1/(delta_t^2)*x(:,i-1)+abar)+C*(1/2/delta_t*x(:,i-1)+abar);
    x(:,i+1) = iK*F_tr(:,i);
    v(:,i+1) = (x(:,i+1)-x(:,i-1))/2/delta_t;
    a(:,i+1) = (x(:,i+1)-2*x(:,i)+x(:,i-1))/(delta_t^2);
end
end