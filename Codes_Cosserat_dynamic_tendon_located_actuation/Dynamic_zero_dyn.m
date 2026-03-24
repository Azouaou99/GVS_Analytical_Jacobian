function [res,Grad] = Dynamic_zero_dyn(var,q0,qd0,Param)
 global  A B C G w
q1=var(1:Param.n*Param.na);
qd1=(q1-q0(1:Param.n*Param.na))/Param.dt;%var(2*Param.n*Param.na+1:3*Param.n*Param.na);
qdd1=(qd1-qd0(1:Param.n*Param.na))/Param.dt;%var(2*Param.n*Param.na+1:3*Param.n*Param.na);
SM1=zeros(Param.n*Param.na,Param.n*Param.na);
SC1=zeros(Param.n*Param.na,Param.n*Param.na);
SFe1=zeros(Param.n*Param.na,1);
SFa1=zeros(Param.n*Param.na,1);
qd1q=zeros(Param.n*Param.na,1);

P_p=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_p=[0;0;0;1];
P_R=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_R=[1 0 0;0 1 0;0 0 1;0 0 0];

P_y=[0 1 0 0];
Q_y=[0; 0; 0;1];

g1=Param.g0;
J1= zeros(6,Param.n*Param.na);
I1= zeros(6,Param.n*Param.na);
J1d= zeros(6,Param.n*Param.na);
I1d= zeros(6,Param.n*Param.na);
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
Pd1=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
phival=Phi(Param.na,Param.n,Y(1),Param.L);
P1=Ad_function(g1)*Param.B*phival;
P1q=zeros(6,Param.n*Param.na*Param.n*Param.na);
N1q=zeros(6,Param.n*Param.na*Param.n*Param.na);
I1q=zeros(6,Param.n*Param.na*Param.n*Param.na);
J1q=zeros(6,Param.n*Param.na*Param.n*Param.na);

P1dq=zeros(6,Param.n*Param.na*Param.n*Param.na);
N1dq=zeros(6,Param.n*Param.na*Param.n*Param.na);
I1dq=zeros(6,Param.n*Param.na*Param.n*Param.na);
J1dq=zeros(6,Param.n*Param.na*Param.n*Param.na);

P1dqd=zeros(6,Param.n*Param.na*Param.n*Param.na);
N1dqd=zeros(6,Param.n*Param.na*Param.n*Param.na);
I1dqd=zeros(6,Param.n*Param.na*Param.n*Param.na);
J1dqd=zeros(6,Param.n*Param.na*Param.n*Param.na);

fM1q=zeros(Param.n*Param.na,Param.n*Param.na);
SM1q=zeros(Param.n*Param.na,Param.n*Param.na);

fC1q=zeros(Param.n*Param.na,Param.n*Param.na);
SC1q=zeros(Param.n*Param.na,Param.n*Param.na);

fC1qd=zeros(Param.n*Param.na,Param.n*Param.na);
SC1qd=zeros(Param.n*Param.na,Param.n*Param.na);

fFe1q=zeros(Param.n*Param.na,Param.n*Param.na);
SFe1q=zeros(Param.n*Param.na,Param.n*Param.na);

H1=Ftendons(Y(1),1,q1,Param);
fFa10=-phival'*Param.B'*H1;
fFe10=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
fM10=J1'*Param.M*J1;
fC10=J1'*(Param.M*J1d -ad_func(J1*qd1)' * Param.M*J1);

fC1qd0=fC1qd;
fM1q0=fM1q;
fC1q0=fC1q;
fFe1q0=fFe1q;
h=Param.dX;
C=zeros(1,2*Param.n*Param.na);
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
    I1d=I1d+(Pd1+Nd1)*Param.dX/2;
    J1d=-Ad_function(Inverse_T(g1))*I1d;
    
    H1=Ftendons(Y(k),k,q1,Param);
    fFa1=-phival'*Param.B'*H1;
    fM1=J1'*Param.M*J1;
    fC1=J1'*(Param.M*J1d -ad_func(J1*qd1)' * Param.M*J1);
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

 for i=1:Param.n*Param.na
        qd1q(i)=1/Param.dt;
        Grad_g1_q=g1*Hat(J1(:,i));
       if k==Param.n_X+1
        C(i)=P_y*Grad_g1_q*Q_y;
       end
        
        Ad_g1_q=[P_R*Grad_g1_q*Q_R zeros(3);(Skew_symmetric(P_p*(g1*Q_p))*(P_R*Grad_g1_q*Q_R))+ (Skew_symmetric(P_p*(Grad_g1_q*Q_p))*(P_R*g1*Q_R)) P_R*Grad_g1_q*Q_R];
        a1=-(Skew_symmetric((P_R*Grad_g1_q*Q_R)'*(P_p*(g1*Q_p))+(P_R*g1*Q_R)'*(P_p*(Grad_g1_q*Q_p)))*(P_R*g1*Q_R)'+Skew_symmetric((P_R*g1*Q_R)'*(P_p*(g1*Q_p)))*(P_R*Grad_g1_q*Q_R)');
        Ad_inv_g1_q=[(P_R*Grad_g1_q*Q_R)' zeros(3);a1 (P_R*Grad_g1_q*Q_R)'];

        %P1q=Ad_inv_g1_q(:,6*(i-1)+1:6*(i))*Param.B*phival;
        N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_g1_q*Param.B*phival;
        I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_inv_g1_q*I1+Ad_function(Inverse_T(g1))*I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i));
  
        N1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_g1_q*ad_func(Param.B*phival*qd1)*J1+ Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+Ad_function(g1)*ad_func(Param.B*phival*qd1q)*J1;
        I1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        J1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=-(Ad_inv_g1_q*I1d+Ad_function(Inverse_T(g1))*I1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)));       
       
        N1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_function(g1)*ad_func(Param.B*phival(:,i))*J1;
        I1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        J1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=-(Ad_function(Inverse_T(g1))*I1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)));

        fM1q(:,i)=2*J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))'*Param.M*J1*qdd1;
        fC1q(:,i)=(J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))'*(Param.M*J1d -ad_func(J1*qd1)' * Param.M*J1)+J1'*(Param.M*J1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))-ad_func( J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))*qd1+J1*qd1q)' * Param.M*J1-ad_func(J1*qd1)' * Param.M*J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))))*qd1;
        fFe1q(:,i)=-( J1'*Param.M*Ad_inv_g1_q*Param.g + J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))'*Param.M*Ad_function(Inverse_T(g1))*Param.g );
        
        fC1qd(:,i)= (J1'*(Param.M*J1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))-ad_func(J1(:,i))' * Param.M*J1))*qd1;
        
        qd1q(i)=0;
 end
    if k==Param.n_X+1
        fM1qn=fM1q;
        fC1qn=fC1q; 
        fFe1qn=fFe1q;
        fC1qdn=fC1qd;
    end
    if and(k~=Param.n_X+1,k~=1)
        SM1q=SM1q+fM1q;
        SC1q=SC1q+fC1q; 
        SFe1q=SFe1q+fFe1q;
        SC1qd=SC1qd+fC1qd;
    end
    P1q=N1q;
    P1dq=N1dq;
    P1dqd=N1dqd;
end
%r1=g1(1:3,4);
%R1=g1(1:3,1:3);

%Fe1=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
M=Param.L/(Param.n_X)*((fM10+fM1n)/2 + SM1);
Co=Param.L/(Param.n_X)*((fC10+fC1n)/2 + SC1);
Mq=Param.L/(Param.n_X)*((fM1q0+fM1qn)/2 + SM1q);
Cq=Param.L/(Param.n_X)*((fC1q0+fC1qn)/2 + SC1q);
Cqd=Param.L/(Param.n_X)*((fC1qd0+fC1qdn)/2 + SC1qd);
Feq=Param.L/(Param.n_X)*((fFe1q0+fFe1qn)/2 + SFe1q);
Fe=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
Hq= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
Keps=Param.Keps;
Deps=Param.Deps;

Grad= M*(Param.dt)^(-2)*eye(Param.n*Param.na)+(Co+Deps)*(Param.dt)^(-1)*eye(Param.n*Param.na)+Mq+Cq+Keps+Feq;
%Grad=Mq+Cq+Keps+Feq;
res=M*qdd1-Hq*Param.Forces_Tendons+Fe+Co*qd1+Keps*q1+Deps*qd1;%2nd derivative of q
K_star=Keps+Feq;
C_star=Cqd+Deps;
M_inv=inv(M);
A=[zeros(Param.n*Param.na), eye(Param.n*Param.na); -M_inv*K_star, -M_inv*C_star];
B=[zeros(Param.n*Param.na,1);M_inv*Hq];
sys = ss(A, B, C, 0);
G = tf(sys);
[V,eigen]= eig(A,"vector");
%[V1,eigen1]= eig(-M_sinv*K_star,"vector");
w=abs(eigen);


