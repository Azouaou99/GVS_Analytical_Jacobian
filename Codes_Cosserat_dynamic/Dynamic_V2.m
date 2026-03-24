function res = Dynamic_V2(var,q0,qd0,i,Param)
q1=var(1:Param.n*Param.na);
qd1=(q1-q0(1:Param.n*Param.na))/Param.dt;%var(2*Param.n*Param.na+1:3*Param.n*Param.na);
qdd1=(qd1-qd0(1:Param.n*Param.na))/Param.dt;%var(2*Param.n*Param.na+1:3*Param.n*Param.na);
SM1=zeros(Param.n*Param.na,Param.n*Param.na);
SC1=zeros(Param.n*Param.na,Param.n*Param.na);
SFe1=zeros(Param.n*Param.na,1);
SFa1=zeros(Param.n*Param.na,1);

g1=Param.g0;
J1= zeros(6,Param.n*Param.na);
I1= zeros(6,Param.n*Param.na);

Jd1= zeros(6,Param.n*Param.na);
Id1= zeros(6,Param.n*Param.na);
% Y1=0:Param.dX:Param.L;
%  for j=1:Param.n_X
%      Jp0=Jp(:,(j-1)*Param.n*Param.na+1:(j)*Param.n*Param.na);
%      Jp(:,(j)*Param.n*Param.na+1:(j+1)*Param.n*Param.na)=Jacobian_subseg(Jp0,Y1(j+1),q,Param);
%      Jp1=Jp(:,(j)*Param.n*Param.na+1:(j+1)*Param.n*Param.na);
%      Jdp0=Jdp(:,(j-1)*Param.n*Param.na+1:(j)*Param.n*Param.na);
%      Jdp(:,(j)*Param.n*Param.na+1:(j+1)*Param.n*Param.na)=dJacobian_subseg(Jp1,Jp0,Jdp0,Y1(j+1),q,qd,Param);
%  end
Y=0:Param.dX:Param.L;
phival=Phi(Param.na,Param.n,Y(1),Param.L);
H1=Ftendons(Y(1),1,q1,Param);
fFa10=-phival'*Param.B'*H1;
fM10=J1'*Param.M*J1;
fC10=J1'*(Param.M*Jd1 -ad_func(J1*qd1)' * Param.M*J1);
fFe10=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
P1=Ad_function(g1)*Param.B*phival;
Pd1=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
h=Param.dX;
for k=2:Param.n_X+1
%     xiaval1=Param.xi_a0+phival*q1;
%     xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
%     phival=Phi(Param.na,Param.n,Y(k),Param.L);
%     g1=g1*expm(Hat(xival1)*(Param.dX));
     %g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));
    x1=Y(k-1)+h/2-sqrt(3)*h/6;
    x2=Y(k-1)+h/2+sqrt(3)*h/6;
    phival1=Phi(Param.na,Param.n,x1,Param.L);
    phival2=Phi(Param.na,Param.n,x2,Param.L);
    xiaval11=Param.xi_a0+phival1*q1;
    xival11=Param.B*xiaval11+Param.B_bar*Param.xi_c;
   
    xiaval12=Param.xi_a0+phival2*q1;
    xival12=Param.B*xiaval12+Param.B_bar*Param.xi_c;

    
    om_h1=h/2*(xival11+xival12)+sqrt(3)*h^2/12* ad_func(xival11)*xival12;   
    g1=g1*expm(Hat(om_h1));
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    
    N1=Ad_function(g1)*Param.B*phival;
    I1=I1+(P1+N1)*Param.dX/2;
    J1=Ad_function(Inverse_T(g1))*I1;
    Nd1=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
    Id1=Id1+(Pd1+Nd1)*Param.dX/2;
    Jd1=-Ad_function(Inverse_T(g1))*Id1;
    
    H1=Ftendons(Y(k),k,q1,Param);
    fFa1=-phival'*Param.B'*H1;
    fM1=J1'*Param.M*J1;
    fC1=J1'*(Param.M*Jd1 -ad_func(J1*qd1)' * Param.M*J1);
    fFe1=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
    %fFa=-phival'*Param.B'*(Ftendonsval);
 
    if k==Param.n_X+1
        fM1n=fM1;
        fC1n=fC1;
        fFe1n=fFe1;
        fFa1n=fFa1;
    end
    if and(k~=Param.n_X+1,k~=1)
        SM1=SM1+fM1;
        SC1=SC1+fC1;
        SFe1=SFe1+fFe1;
        SFa1=SFa1+fFa1;
    end
    %Q1=Q1 + Q1_X*Param.dX;
    %r1=r1 + r1_X*Param.dX;
P1=N1;
Pd1=Nd1;
end
r1=g1(1:3,4);
R1=g1(1:3,1:3);

Mq=Param.L/(Param.n_X)*((fM10+fM1n)/2 + SM1);
Cqqd=Param.L/(Param.n_X)*((fC10+fC1n)/2 + SC1);
Fe=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
Hq= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
Keps=Param.Keps;
Deps=Param.Deps;

res=Mq*qdd1-Hq*Param.Forces_Tendons+Fe+Cqqd*qd1+Keps*q1+Deps*qd1;%2nd derivative of q