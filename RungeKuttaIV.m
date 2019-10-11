function [X,V,A]=RungeKuttaIV(M,K,C,FE,x0,v0,dt,RecordLength)

%=======================================================================
t=(0:dt:(RecordLength-1)*dt);
th=(dt/2:dt:RecordLength*dt);
tf=(dt:dt:RecordLength*dt);
F=Fun(t);
Fh=Fun(th);
Ff=Fun(tf);


X=zeros(size(F,1),RecordLength);
V=zeros(size(F,1),RecordLength);
A=zeros(size(F,1),RecordLength);
X(:,1) = x0;
V(:,1) = v0;
% A(:,1) = inv(M)*(F(:,1)-C*v0-K*x0);
A(:,1)=M\(F(:,1)-C*V(:,1)-K*X(:,1));


for i=1:RecordLength-1
    k11=V(:,i);
    k21=M\(-C*V(:,i)-K*X(:,i)+F(:,i));
    k12=V(:,i)+0.5*dt*k21;
    k22=M\(-C*(V(:,i)+0.5*dt*k21)-K*(X(:,i)+0.5*dt*k11)+Fh(:,i));
    k13=V(:,i)+0.5*dt*k22;
    k23=M\(-C*(V(:,i)+0.5*dt*k22)-K*(X(:,i)+0.5*dt*k12)+Fh(:,i));
    k14=V(:,i)+dt*k23;
    k24=M\(-C*(V(:,i)+dt*k23)-K*(X(:,i)+dt*k13)+Ff(:,i));
    
    X(:,i+1)=X(:,i)+(dt/6)*(k11+2*k12+2*k13+k14);
    V(:,i+1)=V(:,i)+(dt/6)*(k21+2*k22+2*k23+k24);
    A(:,i+1)=M\(F(:,i+1)-C*V(:,i+1)-K*X(:,i+1));
end