function [res,Grad] =Static_V2_grad(q1,Param)
%q1=var(1:Param.n*Param.na);
SFe1=zeros(Param.n*Param.na,1);
SFa1=zeros(Param.n*Param.na,1);

P_p=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_p=[0;0;0;1];
P_R=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_R=[1 0 0;0 1 0;0 0 1;0 0 0];

g1=Param.g0;
J1= zeros(6,Param.n*Param.na);
I1= zeros(6,Param.n*Param.na);

Y=0:Param.dX:Param.L;
phival=Phi(Param.na,Param.n,Y(1),Param.L);
P1=Ad_function(g1)*Param.B*phival;
P1q=zeros(6,Param.n*Param.na*Param.n*Param.na);
N1q=zeros(6,Param.n*Param.na*Param.n*Param.na);
I1q=zeros(6,Param.n*Param.na*Param.n*Param.na);
J1q=zeros(6,Param.n*Param.na*Param.n*Param.na);


fFe1q=zeros(Param.n*Param.na,Param.n*Param.na);
SFe1q=zeros(Param.n*Param.na,Param.n*Param.na);

H1=Ftendons(Y(1),1,q1,Param);
fFa10=-phival'*Param.B'*H1;
fFe10=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
fFe1q0=fFe1q;
h=Param.dX;

for k=2:Param.n_X+1
%      xiaval1=Param.xi_a0+phival*q1;
%      xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
% %     phival=Phi(Param.na,Param.n,Y(k),Param.L);
% %     g1=g1*expm(Hat(xival1)*(Param.dX));
% 
%      g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));

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
    
    H1=Ftendons(Y(k),k,q1,Param);
    fFa1=-phival'*Param.B'*H1;
    fFe1=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
    %fFa=-phival'*Param.B'*(Ftendonsval);
 
    if k==Param.n_X+1
        fFe1n=fFe1;
        fFa1n=fFa1;
    end
    if and(k~=Param.n_X+1,k~=1)
        SFe1=SFe1+fFe1;
        SFa1=SFa1+fFa1;
    end
    %Q1=Q1 + Q1_X*Param.dX;
    %r1=r1 + r1_X*Param.dX;
P1=N1;

 for i=1:Param.n*Param.na
        Grad_g1_q=g1*Hat(J1(:,i));
        Ad_g1_q=[P_R*Grad_g1_q*Q_R zeros(3);(Skew_symmetric(P_p*(g1*Q_p))*(P_R*Grad_g1_q*Q_R))+ (Skew_symmetric(P_p*(Grad_g1_q*Q_p))*(P_R*g1*Q_R)) P_R*Grad_g1_q*Q_R];
        a1=-(Skew_symmetric((P_R*Grad_g1_q*Q_R)'*(P_p*(g1*Q_p))+(P_R*g1*Q_R)'*(P_p*(Grad_g1_q*Q_p)))*(P_R*g1*Q_R)'+Skew_symmetric((P_R*g1*Q_R)'*(P_p*(g1*Q_p)))*(P_R*Grad_g1_q*Q_R)');
        Ad_inv_g1_q=[(P_R*Grad_g1_q*Q_R)' zeros(3);a1 (P_R*Grad_g1_q*Q_R)'];

        %P1q=Ad_inv_g1_q(:,6*(i-1)+1:6*(i))*Param.B*phival;
        N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_g1_q*Param.B*phival;
        I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_inv_g1_q*I1+Ad_function(Inverse_T(g1))*I1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i));
        fFe1q(:,i)=-( J1'*Param.M*Ad_inv_g1_q*Param.g + J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))'*Param.M*Ad_function(Inverse_T(g1))*Param.g );
 end
    if k==Param.n_X+1 
        fFe1qn=fFe1q;
    else
        SFe1q=SFe1q+fFe1q;
    end
    P1q=N1q;
end
r1=g1(1:3,4);
R1=g1(1:3,1:3);

%Fe1=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);

Feq=Param.dX.*((fFe1q0+fFe1qn)/2 + SFe1q);
Fe=Param.dX.*((fFe10+fFe1n)/2 + SFe1);

Hq= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
Keps=Param.Keps;

Grad= Keps+Feq;
%Fep=[0;0.01;0;0;0;0];
%Grad=Mq+Cq+Keps+Feq;
%J1'
res=-Hq*Param.Forces_Tendons+Fe-J1'*Param.Ftip+Keps*q1;%2nd derivative of q
