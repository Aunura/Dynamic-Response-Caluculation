function [X,V,A]=Houbolt(M,K,C,F,X0,V0,dt,RecordLength)
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
%       The second and third step will be calculated by CentralDifference
%      
%================================================================
[XI,VI,AI]=CentralDifferenceM(M,K,C,F(:,1:3),dt,X0,V0,3);

X=zeros(size(F,1),RecordLength);
V=zeros(size(F,1),RecordLength);
A=zeros(size(F,1),RecordLength);

for i=1:3
   X(:,i)=XI(:,i); 
   V(:,i)=VI(:,i); 
   A(:,i)=AI(:,i); 
end


a1=2/(dt^2);a2=12/(6*dt);
b1=5/(dt^2);b2=3/dt;b3=4/(dt^2);b4=3/(2*dt);b5=1/(dt^2);b6=1/(3*dt);
c1=1/(dt^2);c2=1/(6*dt);

ks=a1*M+a2*C+K;


for i=3:RecordLength-1
    Fs=F(:,i+1)+(b1*M+b2*C)*X(:,i)-(b3*M+b4*C)*X(:,i-1)+(b5*M+b6*C)*X(:,i-2);
    X(:,i+1)=ks\Fs;
    A(:,i+1)=c1*(2*X(:,i+1)-5*X(:,i)+4*X(:,i-1)-X(:,i-2));
    V(:,i+1)=c2*(11*X(:,i+1)-18*X(:,i)+9*X(:,i-1)-2*X(:,i-2));
end

end


function [x,v,a]=CentralDifferenceM(m,k,c,f,dt,x0,v0,recordlenght)
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
a(:,1) = inv(m)*(f(:,1)-c*v0-k*x0);
a0 = a(:,1);

K_tr = m/(delta_t^2)+c/2/delta_t;
iK=inv(K_tr);

%=============解决起步问题（前二个时间点的求解）===================
x_neg1 = x0-v0*delta_t+(a0/2)*delta_t^2;
F_tr(:,1) = f(:,1)-(m/(delta_t^2)-c/2/delta_t)*x_neg1-(k-2*m/(delta_t^2))*x(:,1);
x(:,2) = iK*F_tr(:,1);
v(:,2) = (x(:,2)-x(:,1))/2/delta_t;
a(:,2) = (x(:,2)-2*x(:,1)+x_neg1)/(delta_t^2);
%===================end========================================

for i=2:recordlenght-1
    F_tr(:,i) = f(:,i)-(m/(delta_t^2)-c/2/delta_t)*x(:,i-1)-(k-2*m/(delta_t^2))*x(:,i);
%     F_tr(:,i) = F(:,i)-K*(x(:,i)+abar)+M*(2*x(:,i)/(delta_t^2)-1/(delta_t^2)*x(:,i-1)+abar)+C*(1/2/delta_t*x(:,i-1)+abar);
    x(:,i+1) = iK*F_tr(:,i);
    v(:,i+1) = (x(:,i+1)-x(:,i-1))/2/delta_t;
    a(:,i+1) = (x(:,i+1)-2*x(:,i)+x(:,i-1))/(delta_t^2);
end
end