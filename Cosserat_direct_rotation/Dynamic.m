function [res,Grad] = Dynamic(var,q0,qd0,qa0,eta0,i,Param)
 global  A B C G w
qa1=zeros(6,1);
qa1(3)=var(1);
q1=var(2:Param.n*Param.na+1);
q=[qa1;q1];
qd1=(q1-q0(1:Param.n*Param.na))/Param.dt;%var(Param.n*Param.na+1:3*Param.n*Param.na);
qdd1=(qd1-qd0(1:Param.n*Param.na))/Param.dt;%var(Param.n*Param.na+1:3*Param.n*Param.na);

eta01=(qa1-qa0(1:6))/Param.dt;
eta01d=(eta01-eta0(1:6))/Param.dt;
qd=[eta01;qd1];
qdd=[eta01d;qdd1];

qd1q=zeros(6+Param.n*Param.na,1);
qld1q=zeros(Param.n*Param.na,1);
SM1=zeros(Param.n*Param.na+6,Param.n*Param.na+6);
SC1=zeros(Param.n*Param.na+6,Param.n*Param.na+6);
SFe1=zeros(Param.n*Param.na+6,1);
SFa1=zeros(Param.n*Param.na,1);

R01=[cos(qa1(3)) -sin(qa1(3)) 0; sin(qa1(3)) cos(qa1(3)) 0; 0 0 1];
g10=Param.g0;
g10(1:3,1:3)=R01;


J1= zeros(6,Param.n*Param.na);
J10= zeros(6,6);

I1= zeros(6,Param.n*Param.na);

J1d= zeros(6,Param.n*Param.na);
I1d= zeros(6,Param.n*Param.na);
I1d0= zeros(6,6);
J1d0= zeros(6,6);
P_p=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_p=[0;0;0;1];
P_R=[1 0 0 0;0 1 0 0;0 0 1 0];
Q_R=[1 0 0;0 1 0;0 0 1;0 0 0];

P_y=[0 1 0 0];
Q_y=[0; 0; 0;1];

g1=g10;
Y=0:Param.dX:Param.L;
% Grad_g1_q=zeros(4,4*(6+Param.n*Param.na));
% 
% Grad_g1_qq=zeros(4,4*(6+Param.n*Param.na)*(6+Param.n*Param.na));
% Grad_g2_qq=zeros(4,4*(6+Param.n*Param.na)*(6+Param.n*Param.na));
phival=Phi(Param.na,Param.n,Y(1),Param.L);
P1=Ad_function(g1)*Param.B*phival;
%P2=Ad_function(g2)*Param.B*phival;
P1q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%P2q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
N1q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%N2q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
I1q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%I2q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
J1q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%J2q=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
J1q0=zeros(6,6*(Param.n*Param.na+6));
%J2q0=zeros(6,6*(Param.n*Param.na+6));

P1d0=zeros(6);
%P2d0=zeros(6);
P1d=zeros(6,Param.n*Param.na);
%P2d=zeros(6,Param.n*Param.na);

P1dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%P2dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
P1dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%P2dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));

N1dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%N2dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
N1dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%N2dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));

I1dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%I2dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
I1dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%I2dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));

J1dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%J2dq=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
J1dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));
%J2dqd=zeros(6,Param.n*Param.na*(Param.n*Param.na+6));


P1dq0=zeros(6,6*(Param.n*Param.na+6));
%P2dq0=zeros(6,6*(Param.n*Param.na+6));
N1dq0=zeros(6,6*(Param.n*Param.na+6));
%N2dq0=zeros(6,6*(Param.n*Param.na+6));
I1dq0=zeros(6,6*(Param.n*Param.na+6));
%I2dq0=zeros(6,6*(Param.n*Param.na+6));
J1dq0=zeros(6,6*(Param.n*Param.na+6));
%J2dq0=zeros(6,6*(Param.n*Param.na+6));

P1dqd0=zeros(6,6*(Param.n*Param.na+6));
%P2dqd0=zeros(6,6*(Param.n*Param.na+6));
N1dqd0=zeros(6,6*(Param.n*Param.na+6));
%N2dqd0=zeros(6,6*(Param.n*Param.na+6));
I1dqd0=zeros(6,6*(Param.n*Param.na+6));
%I2dqd0=zeros(6,6*(Param.n*Param.na+6));
J1dqd0=zeros(6,6*(Param.n*Param.na+6));
%J2dqd0=zeros(6,6*(Param.n*Param.na+6));

fM1q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%fM2q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
SM1q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%SM2q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);

fC1q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%fC2q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
SC1q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%SC2q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);

fC1qd=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%fC2qd=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
SC1qd=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%SC2qd=zeros(6+Param.n*Param.na,6+Param.n*Param.na);

fFe1q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
SFe1q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%fFe2q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);
%SFe2q=zeros(6+Param.n*Param.na,6+Param.n*Param.na);

H1=Ftendons(Y(1),1,q1,Param);
%H2=Ftendons(Y(1),1,q2,Param);
fFa10=-phival'*Param.B'*H1;
%fFa20=-phival'*Param.B'*H2;
fFe10=[-J10'*Param.M*Ad_function(Inverse_T(g1))*Param.g;-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g];
%fFe20=[-J20'*Param.M*Ad_function(Inverse_T(g2))*Param.g;-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g];

fM10=[J10';J1']*Param.M*[J10 J1];
%fM20=[J20';J2']*Param.M*[J20 J2];

fC10=[J10';J1']*(Param.M*[J1d0 J1d] -ad_func([J10 J1]*[eta01;qd1])' * Param.M*[J10 J1]);
%fC20=[J20';J2']*(Param.M*[J2d0 J2d] -ad_func([J20 J2]*[eta02;qd2])' * Param.M*[J20 J2]);

fM1q0=fM1q;
%fM2q0=fM2q;
fC1q0=fC1q;
%fC2q0=fC2q;
fC1qd0=fC1qd;
%fC2qd0=fC2qd;
fFe1q0=fFe1q;
%fFe2q0=fFe2q;
h=Param.dX;
C=zeros(1,2*Param.n*Param.na+12);
for k=2:Param.n_X+1
%     xiaval1=Param.xi_a0+phival*q1;
%     xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
%     xiaval2=Param.xi_a0+phival*q2;
%     xival2=Param.B*xiaval2+Param.B_bar*Param.xi_c;
%     phival=Phi(Param.na,Param.n,Y(k),Param.L);
%     g1=g1*expm(Hat(xival1)*(Param.dX));
%     g2=g2*expm(Hat(xival2)*(Param.dX));
     %g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));
     %g2=g2*expm2(Hat(xival2)*(Param.dX),xival2*(Param.dX));
    x1=Y(k-1)+h/2-sqrt(3)*h/6;
   x2=Y(k-1)+h/2+sqrt(3)*h/6;
    phival1=Phi(Param.na,Param.n,x1,Param.L);
    phival2=Phi(Param.na,Param.n,x2,Param.L);
    xiaval11=Param.xi_a0+phival1*q1;
    xival11=Param.B*xiaval11+Param.B_bar*Param.xi_c;
   % xiaval21=Param.xi_a0+phival1*q2;
    %xival21=Param.B*xiaval21+Param.B_bar*Param.xi_c;
   
    xiaval12=Param.xi_a0+phival2*q1;
    xival12=Param.B*xiaval12+Param.B_bar*Param.xi_c;
    %xiaval22=Param.xi_a0+phival2*q2;
    %xival22=Param.B*xiaval22+Param.B_bar*Param.xi_c;
    
    om_h1=h/2*(xival11+xival12)+sqrt(3)*h^2/12* ad_func(xival11)*xival12;
  %  om_h2=h/2*(xival21+xival22)+sqrt(3)*h^2/12* ad_func(xival21)*xival22;    
    g1=g1*expm(Hat(om_h1));
   % g2=g2*expm(Hat(om_h2));
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    
    N1=Ad_function(g1)*Param.B*phival;
    %N2=Ad_function(g2)*Param.B*phival;
    I1=I1+(P1+N1)*Param.dX/2;
    %I2=I2+(P2+N2)*Param.dX/2;
    J1=Ad_function(Inverse_T(g1))*I1;
    %J2=Ad_function(Inverse_T(g2))*I2;

    J10=Ad_function(Inverse_T(g1))*Ad_function(g10);
    %J20=Ad_function(Inverse_T(g2))*Ad_function(g20);

    N1d=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
    %N2d=Ad_function(g2)*ad_func(Param.B*phival*qd2)*J2;
    I1d=I1d+(P1d+N1d)*Param.dX/2;
    %I2d=I2d+(P2d+N2d)*Param.dX/2;
    J1d=-Ad_function(Inverse_T(g1))*I1d;
    %J2d=-Ad_function(Inverse_T(g2))*I2d;

    N1d0=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J10;
    %N2d0=Ad_function(g2)*ad_func(Param.B*phival*qd2)*J20;
    I1d0=I1d0+(P1d0+N1d0)*Param.dX/2;
    %I2d0=I2d0+(P2d0+N2d0)*Param.dX/2;
    J1d0=-Ad_function(Inverse_T(g1))*I1d0;
    %J2d0=-Ad_function(Inverse_T(g2))*I2d0;



    H1=Ftendons(Y(k),k,q1,Param);
    %H2=Ftendons(Y(k),k,q2,Param);
    fFa1=-phival'*Param.B'*H1;
   %fFa2=-phival'*Param.B'*H2;
    fFe1=[-J10'*Param.M*Ad_function(Inverse_T(g1))*Param.g;-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g];
    %fFe2=[-J20'*Param.M*Ad_function(Inverse_T(g2))*Param.g;-J2'*Param.M*Ad_function(Inverse_T(g2))*Param.g];
    
    fM1=[J10';J1']*Param.M*[J10 J1];
    %fM2=[J20';J2']*Param.M*[J20 J2];

    fC1=[J10';J1']*(Param.M*[J1d0 J1d] -ad_func([J10 J1]*[eta01;qd1])' * Param.M*[J10 J1]);
    %fC2=[J20';J2']*(Param.M*[J2d0 J2d] -ad_func([J20 J2]*[eta02;qd2])' * Param.M*[J20 J2]);
    %fFa=-phival'*Param.B'*(Ftendonsval);
    if k==Param.n_X+1
        fFe1n=fFe1;
     %   fFe2n=fFe2;
        fFa1n=fFa1;
     %   fFa2n=fFa2;
        
        fM1n=fM1;
     %   fM2n=fM2;
        fC1n=fC1;
     %   fC2n=fC2;
    end
    if and(k~=Param.n_X+1,k~=1)
        SFe1=SFe1+fFe1;
      %  SFe2=SFe2+fFe2;
        SFa1=SFa1+fFa1;
      %  SFa2=SFa2+fFa2;
        SM1=SM1+fM1;
      %  SM2=SM2+fM2;
        SC1=SC1+fC1;
      %  SC2=SC2+fC2;
    end
    P1=N1;
    %P2=N2;
    P1d=N1d;
    %P2d=N2d;
    P1d0=N1d0;
    %P2d0=N2d0;

    for i=1:6+Param.n*Param.na
        qd1q(i)=Param.dt^(-1);
     %   qd2q(i)=Param.dt^(-1);
        if i>6
        qld1q(i-6)=Param.dt^(-1);
      %  qld2q(i-6)=Param.dt^(-1);
        end
       
        if i >= 1 && i <= 6
        Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J10(:,i));
       % Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))=g2*Hat(J20(:,i));
        else
        Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i-6));
        %Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))=g2*Hat(J2(:,i-6));
        end
       if k==Param.n_X+1
        C(i)=P_y*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_y;
       end
        
        %Grad_g1_q(:,4*(i-1)+1:4*(i))=g1*Hat(J1(:,i));
        %Grad_g2_q(:,4*(Param.n*Param.na+i-1)+1:4*(Param.n*Param.na+i))=g2*Hat(J2(:,i));

        Ad_g1_q=[P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R zeros(3);(Skew_symmetric(P_p*(g1*Q_p))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p))*(P_R*g1*Q_R)) P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R];
        %Ad_g2_q=[P_R*Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_R zeros(3);(Skew_symmetric(P_p*(g2*Q_p))*(P_R*Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_R))+ (Skew_symmetric(P_p*(Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_p))*(P_R*g2*Q_R)) P_R*Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_R];
        a1=-(Skew_symmetric((P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'*(P_p*(g1*Q_p))+(P_R*g1*Q_R)'*(P_p*(Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_p)))*(P_R*g1*Q_R)'+Skew_symmetric((P_R*g1*Q_R)'*(P_p*(g1*Q_p)))*(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)');
        %a2=-(Skew_symmetric((P_R*Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_R)'*(P_p*(g2*Q_p))+(P_R*g2*Q_R)'*(P_p*(Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_p)))*(P_R*g2*Q_R)'+Skew_symmetric((P_R*g2*Q_R)'*(P_p*(g2*Q_p)))*(P_R*Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_R)');
        Ad_inv_g1_q=[(P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)' zeros(3);a1 (P_R*Grad_g1_q(:,4*(i-1)+1:4*(i))*Q_R)'];
        %Ad_inv_g2_q=[(P_R*Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_R)' zeros(3);a2 (P_R*Grad_g2_q(:,4*(6+Param.n*Param.na+i-1)+1:4*(6+Param.n*Param.na+i))*Q_R)'];
 % 
        N1q(:,(Param.n*Param.na)*(i-1)+1:(Param.n*Param.na)*(i))=Ad_g1_q*Param.B*phival;
        I1q(:,(Param.n*Param.na)*(i-1)+1:(Param.n*Param.na)*(i))=I1q(:,(Param.n*Param.na)*(i-1)+1:(Param.n*Param.na)*(i))+(P1q(:,(Param.n*Param.na)*(i-1)+1:(Param.n*Param.na)*(i))+N1q(:,(Param.n*Param.na)*(i-1)+1:(Param.n*Param.na)*(i)))*Param.dX/2;
        J1q(:,(Param.n*Param.na)*(i-1)+1:(Param.n*Param.na)*(i))=Ad_inv_g1_q*I1+Ad_function(Inverse_T(g1))*I1q(:,(Param.n*Param.na)*(i-1)+1:(Param.n*Param.na)*(i));
      
        %N2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=Ad_g2_q*Param.B*phival;
        %I2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=I2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))+(P2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))+N2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i)))*Param.dX/2;
        %J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=Ad_inv_g2_q*I2 + Ad_function(Inverse_T(g2))*I2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i));
       
        if i>6
        N1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_function(g1)*ad_func(Param.B*phival(:,i-6))*J1;
        I1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        J1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=-(Ad_function(Inverse_T(g1))*I1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)));

        %N2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=Ad_function(g2)*ad_func(Param.B*phival(:,i-6))*J2;
        %I2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=I2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))+(P2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))+N2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i)))*Param.dX/2;
        %J2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=-(Ad_function(Inverse_T(g2))*I2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i)));
        else
        J1dqd(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=zeros(6,Param.n*Param.na);
        %J2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=zeros(6,Param.n*Param.na);
        end

        %Rp_qp1_qp1=[-cos(qa(3)) sin(qa(3)) 0; -sin(qa(3)) -cos(qa(3)) 0; 0 0 0];
        if i==3
        grad_g01_qa=[-sin(qa1(3)) -cos(qa1(3)) 0 0; cos(qa1(3)) -sin(qa1(3)) 0 0; 0 0 0 0; 0 0 0 0];
        %grad_g02_qa=[-sin(qa2(3)) -cos(qa2(3)) 0 0; cos(qa2(3)) -sin(qa2(3)) 0 0; 0 0 0 0; 0 0 0 0];
        % a10=-(Skew_symmetric((P_R*grad_g01_qa*Q_R)'*(P_p*(g10*Q_p))+(P_R*g10*Q_R)'*(P_p*(grad_g01_qa*Q_p)))*(P_R*g10*Q_R)'+Skew_symmetric((P_R*g10*Q_R)'*(P_p*(g10*Q_p)))*(P_R*grad_g01_qa*Q_R)');
        % a20=-(Skew_symmetric((P_R*grad_g02_qa*Q_R)'*(P_p*(g20*Q_p))+(P_R*g20*Q_R)'*(P_p*(grad_g02_qa*Q_p)))*(P_R*g20*Q_R)'+Skew_symmetric((P_R*g20*Q_R)'*(P_p*(g20*Q_p)))*(P_R*grad_g02_qa*Q_R)');
        % Ad_inv_g10_qa=[(P_R*grad_g01_qa*Q_R)' zeros(3);a10 (P_R*grad_g01_qa*Q_R)'];
        % Ad_inv_g20_qa=[(P_R*grad_g02_qa*Q_R)' zeros(3);a20 (P_R*grad_g02_qa*Q_R)'];
        Ad_g10_qa=[P_R* grad_g01_qa*Q_R zeros(3);(Skew_symmetric(P_p*(g10*Q_p))*(P_R* grad_g01_qa*Q_R))+ (Skew_symmetric(P_p*( grad_g01_qa*Q_p))*(P_R*g10*Q_R)) P_R* grad_g01_qa*Q_R];
        %Ad_g20_qa=[P_R*grad_g02_qa*Q_R zeros(3);(Skew_symmetric(P_p*(g20*Q_p))*(P_R*grad_g02_qa*Q_R))+ (Skew_symmetric(P_p*(grad_g02_qa*Q_p))*(P_R*g20*Q_R)) P_R*grad_g02_qa*Q_R];
        J1q0(:,6*(i-1)+1:6*(i))=Ad_inv_g1_q*Ad_function(g10) + Ad_function(Inverse_T(g1))*Ad_g10_qa;
        %J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=Ad_inv_g2_q*Ad_function(g20) + Ad_function(Inverse_T(g2))*Ad_g20_qa;
        else
        J1q0(:,6*(i-1)+1:6*(i))=Ad_inv_g1_q*Ad_function(g10); 
        %J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=Ad_inv_g2_q*Ad_function(g20);
        end
        N1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=Ad_g1_q*ad_func(Param.B*phival*qd1)*J1+ Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+Ad_function(g1)*ad_func(Param.B*phival*qld1q)*J1;
        I1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=I1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+(P1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))+N1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)))*Param.dX/2;
        J1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))=-(Ad_inv_g1_q*I1d+Ad_function(Inverse_T(g1))*I1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i)));       
       
        
        %N2dq(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=Ad_g2_q*ad_func(Param.B*phival*qd2)*J2 + Ad_function(g2)*ad_func(Param.B*phival*qd2)*J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))+Ad_function(g2)*ad_func(Param.B*phival*qld2q)*J2;
        %I2dq(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=I2dq(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))+(P2dq(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))+N2dq(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i)))*Param.dX/2;
        %J2dq(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))=-(Ad_inv_g2_q*I2d + Ad_function(Inverse_T(g2))*I2dq(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i)));
       
        if i==3
        N1dq0(:,6*(i-1)+1:6*(i))=Ad_g1_q*ad_func(Param.B*phival*qd1)*J10 + Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1q0(:,6*(i-1)+1:6*(i));
        I1dq0(:,6*(i-1)+1:6*(i))=I1dq0(:,6*(i-1)+1:6*(i))+(P1dq0(:,6*(i-1)+1:6*(i))+N1dq0(:,6*(i-1)+1:6*(i)))*Param.dX/2;
        J1dq0(:,6*(i-1)+1:6*(i))=-(Ad_inv_g1_q*I1d0+Ad_function(Inverse_T(g1))*I1dq0(:,6*(i-1)+1:6*(i)));       
       
      
        J1dqd0(:,6*(i-1)+1:6*(i))=zeros(6);

        %N2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=Ad_g2_q*ad_func(Param.B*phival*qd2)*J20 + Ad_function(g2)*ad_func(Param.B*phival*qd2)*J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i));
        %I2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=I2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))+(P2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))+N2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)))*Param.dX/2;
        %J2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=-(Ad_inv_g2_q*I2d0 + Ad_function(Inverse_T(g2))*I2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)));
        

        %J2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=zeros(6);

        else
        N1dq0(:,6*(i-1)+1:6*(i))=Ad_g1_q*ad_func(Param.B*phival*qd1)*J10+ Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1q0(:,6*(i-1)+1:6*(i))+Ad_function(g1)*ad_func(Param.B*phival*qld1q)*J10;
        I1dq0(:,6*(i-1)+1:6*(i))=I1dq0(:,6*(i-1)+1:6*(i))+(P1dq0(:,6*(i-1)+1:6*(i))+N1dq0(:,6*(i-1)+1:6*(i)))*Param.dX/2;
        J1dq0(:,6*(i-1)+1:6*(i))=-(Ad_inv_g1_q*I1d0+Ad_function(Inverse_T(g1))*I1dq0(:,6*(i-1)+1:6*(i)));       
       
        N1dqd0(:,6*(i-1)+1:6*(i))=Ad_function(g1)*ad_func(Param.B*phival)*J10;
        I1dqd0(:,6*(i-1)+1:6*(i))=I1dqd0(:,6*(i-1)+1:6*(i))+(P1dqd0(:,6*(i-1)+1:6*(i))+N1dqd0(:,6*(i-1)+1:6*(i)))*Param.dX/2;
        J1dqd0(:,6*(i-1)+1:6*(i))=-(Ad_function(Inverse_T(g1))*I1dqd0(:,6*(i-1)+1:6*(i)));

        %N2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=Ad_g2_q*ad_func(Param.B*phival*qd2)*J20 + Ad_function(g2)*ad_func(Param.B*phival*qd2)*J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))+Ad_function(g2)*ad_func(Param.B*phival*qld2q)*J20;
        %I2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=I2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))+(P2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))+N2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)))*Param.dX/2;
        %J2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=-(Ad_inv_g2_q*I2d0 + Ad_function(Inverse_T(g2))*I2dq0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)));
        
        %N2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=Ad_function(g2)*ad_func(Param.B*phival)*J20;
        %I2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=I2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))+(P2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))+N2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)))*Param.dX/2;
        %J2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))=-(Ad_function(Inverse_T(g2))*I2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)));
        end
        
        fM1q(:,i)=2*[J1q0(:,6*(i-1)+1:6*(i))';J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))']*Param.M*[J10 J1]*[eta01;qdd1];
        %fM2q(:,i)=2*[J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i))';J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))']*Param.M*[J20 J2]*[eta02;qdd2];                                                                                                                                                            
        fC1q(:,i)=([J1q0(:,6*(i-1)+1:6*(i))';J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))']*(Param.M*[J1d0 J1d] -ad_func([J10 J1]*[eta01;qd1])'* Param.M*[J10 J1])+[J10';J1']*(Param.M*[J1dq0(:,6*(i-1)+1:6*(i)) J1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))]-ad_func( [J1q0(:,6*(i-1)+1:6*(i)) J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))]*[eta01; qd1]+[J10 J1]*qd1q)' * Param.M*[J10 J1]-ad_func([J10 J1]*[eta01;qd1])' * Param.M*[J1q0(:,6*(i-1)+1:6*(i)) J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))]))*[eta01;qd1];
        %fC2q(:,i)=([J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)) J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))]'*(Param.M*[J2d0 J2d] -ad_func([J20 J2]*[eta02;qd2])' * Param.M*[J20 J2])+[J20 J2]'*(Param.M*[J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)) J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))]-ad_func( [J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)) J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))]*[eta02;qd2] + [J20 J2]*qd2q)' * Param.M*[J20 J2]-ad_func([J20 J2]*[eta02;qd2])' * Param.M*[J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)) J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))]))*[eta02;qd2];  
        
        fC1qd(:,i)=([J10';J1']*(Param.M*[J1dqd0(:,6*(i-1)+1:6*(i)) J1dq(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))] -ad_func([J10 J1])' * Param.M*[J10 J1]))*[eta01;qd1];
        %fC2qd(:,i)=([J20';J2']*(Param.M*[J2dqd0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)) J2dqd(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))] -ad_func([J20 J2])' * Param.M*[J20 J2]))*[eta02;qd2];
        
        fFe1q(:,i)=-( [J10 J1]'*Param.M*Ad_inv_g1_q*Param.g + [J1q0(:,6*(i-1)+1:6*(i))';J1q(:,Param.n*Param.na*(i-1)+1:Param.n*Param.na*(i))']*Param.M*Ad_function(Inverse_T(g1))*Param.g );
        %fFe2q(:,i)=-( [J20 J2]'*Param.M*Ad_inv_g2_q*Param.g + [J2q0(:,6*(6+Param.n*Param.na+i-1)+1:6*(6+Param.n*Param.na+i)) J2q(:,Param.n*Param.na*(6+Param.n*Param.na+i-1)+1:Param.n*Param.na*(6+Param.n*Param.na+i))]'*Param.M*Ad_function(Inverse_T(g2))*Param.g );
      
        
       qd1q(i)=0;
       %qd2q(i)=0;
       if i>6
       qld1q(i-6)=0;
       %qld2q(i-6)=0;
       end
    end
    fC1qd=[J10';J1']*(Param.M*[J1d0 J1d] -ad_func([J10 J1]*[eta01;qd1])' * Param.M*[J10 J1])+fC1qd;
    %fC2qd=[J20';J2']*(Param.M*[J2d0 J2d] -ad_func([J20 J2]*[eta02;qd2])' * Param.M*[J20 J2])+fC2qd;
        if k==Param.n_X+1
        fM1qn=fM1q;
     %   fM2qn=fM2q;
        fC1qn=fC1q;
     %   fC2qn=fC2q;
        fFe1qn=fFe1q;
     %   fFe2qn=fFe2q;
        fC1qdn=fC1qd;
     %   fC2qdn=fC2qd;
    end
    if and(k~=Param.n_X+1,k~=1)
        SM1q=SM1q+fM1q;
     %   SM2q=SM2q+fM2q;
        SC1q=SC1q+fC1q;
     %   SC2q=SC2q+fC2q; 
        SC1qd=SC1qd+fC1qd;
     %   SC2qd=SC2qd+fC2qd;
        SFe1q=SFe1q+fFe1q;
     %   SFe2q=SFe2q+fFe2q;
    end

    P1q=N1q;
    %P2q=N2q;
    P1dq=N1dq;
    %P2dq=N2dq;
    P1dqd=N1dqd;
    %P2dqd=N2dqd;
    P1dqd0=N1dqd0;
    %P2dqd0=N2dqd0;
    P1dq0=N1dq0;
    %P2dq0=N2dq0;

end
Feq=Param.L/(Param.n_X)*((fFe1q0+fFe1qn)/2 + SFe1q);
%Fe2q=Param.L/(Param.n_X)*((fFe2q0+fFe2qn)/2 + SFe2q);
%Feq=[Fe1q, zeros(6+Param.n*Param.na,6+Param.n*Param.na+3);zeros(6+Param.n*Param.na), Fe2q zeros(6+Param.n*Param.na,3);zeros(3,12+Param.n*Param.na+3)];

Fe=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
%Fe2=Param.L/(Param.n_X)*((fFe20+fFe2n)/2 + SFe2);

M=Param.L/(Param.n_X)*((fM10+fM1n)/2 + SM1);
%M2=Param.L/(Param.n_X)*((fM20+fM2n)/2 + SM2);
%M=[M1, zeros(Param.n*Param.na+6,3+Param.n*Param.na+6);zeros(Param.n*Param.na+6,Param.n*Param.na+6), M2 , zeros(Param.n*Param.na+6,3); zeros(3,Param.n*Param.na+6) Param.Mp];

Mq=Param.L/(Param.n_X)*((fM1q0+fM1qn)/2 + SM1q);
%M2q=Param.L/(Param.n_X)*((fM2q0+fM2qn)/2 + SM2q);
%Mq=[M1q, zeros(6+Param.n*Param.na,6+Param.n*Param.na+3);zeros(6+Param.n*Param.na), M2q, zeros(6+Param.n*Param.na,3);zeros(3,12+Param.n*Param.na+3)];

Co=Param.L/(Param.n_X)*((fC10+fC1n)/2 + SC1);
%C2=Param.L/(Param.n_X)*((fC20+fC2n)/2 + SC2);
%C=[C1, zeros(Param.n*Param.na+6,3+Param.n*Param.na+6);zeros(Param.n*Param.na+6,Param.n*Param.na+6), C2 , zeros(Param.n*Param.na+6,3); zeros(3,Param.n*Param.na+3+12)];


Cq=Param.L/(Param.n_X)*((fC1q0+fC1qn)/2 + SC1q);
%C2q=Param.L/(Param.n_X)*((fC2q0+fC2qn)/2 + SC2q);
%Cq=[C1q, zeros(6+Param.n*Param.na,6+Param.n*Param.na+3);zeros(6+Param.n*Param.na), C2q zeros(6+Param.n*Param.na,3);zeros(3,12+Param.n*Param.na+3)];

Cqd=Param.L/(Param.n_X)*((fC1qd0+fC1qdn)/2 + SC1qd);
%C2qd=Param.L/(Param.n_X)*((fC2qd0+fC2qdn)/2 + SC2qd);
%Cqd=[C1qd, zeros(6+Param.n*Param.na,6+Param.n*Param.na+3);zeros(6+Param.n*Param.na), C2qd zeros(6+Param.n*Param.na,3);zeros(3,12+Param.n*Param.na+3)];

%gp=[Rp,rp;0 0 0 1];
%Fep=-Param.Mp*Ad_function(Inverse_T(gp))*Param.g;
%Fep=-Param.Mp*Param.g;
%Fep=[0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0]*Fep;
%Fep=[0;-Param.mp*Param.g(4:5)];

%Fep=Param.Ftip+[0;-Param.mp*Param.g(4);0];
%Fe=[Fe1;Fe2;Fep];

Hq= Param.dX.*((fFa10+fFa1n)./2 + SFa1);
%H2q= Param.dX.*((fFa20+fFa2n)./2 + SFa2);
%Hq=[H1q, zeros(Param.n*Param.na,2);zeros(Param.n*Param.na,2), H2q; zeros(3,4)];

Keps=[zeros(6,Param.n*Param.na+6);zeros(Param.n*Param.na,6) Param.Keps];%, zeros(Param.n*Param.na,Param.n*Param.na+3+6) ;zeros(6,Param.n*Param.na+6+3);zeros(Param.n*Param.na,Param.n*Param.na+12),Param.Keps zeros(Param.n*Param.na,3);zeros(3,Param.n*Param.na+3+12)];
Deps=[zeros(6,Param.n*Param.na+6);zeros(Param.n*Param.na,6) Param.Deps];%Deps=[zeros(6,Param.n*Param.na+6+3);zeros(Param.n*Param.na,6) Param.Deps, zeros(Param.n*Param.na,Param.n*Param.na+3+6) ;zeros(6,Param.n*Param.na+6+3);zeros(Param.n*Param.na,Param.n*Param.na+12),Param.Deps zeros(Param.n*Param.na,3);zeros(3,Param.n*Param.na+3+12)];
%Deps=[Param.Deps, zeros(Param.n*Param.na,Param.n*Param.na+3);zeros(Param.n*Param.na), Param.Deps zeros(Param.n*Param.na,3);zeros(3,Param.n*Param.na+3)];
%Rc1=[cos(Param.c1) -sin(Param.c1) 0; sin(Param.c1) cos(Param.c1) 0; 0 0 1];
%Rc2=[cos(Param.c2) -sin(Param.c2) 0; sin(Param.c2) cos(Param.c2) 0; 0 0 1];
%Psi=[X1*(R1'*(r2-r1));X2*(R2'*(r2-r1));sqrt((r1-r2)'*(r1-r2))-Param.y2];


%Grad= [M*(Param.dt)^(-2)*eye(Param.n*Param.na+3)+(C+Deps)*(Param.dt)^(-1)*eye(Param.n*Param.na+3)+Mq+Cq+Keps + Feq+ Psi_q_q_transpose_lambda Psi_q'; Psi_q zeros(4,4)];
%Grad=0;
Grad= [M*(Param.dt)^(-2)*eye(6+Param.n*Param.na)+(Co+Deps)*(Param.dt)^(-1)*eye(6+Param.n*Param.na)+Mq+Cq+Keps + Feq];

Grad=[Grad(3,:);Grad(7:end,:)];
Grad=[Grad(:,3) Grad(:,7:end)];

Act=0*Hq*Param.Forces_Tendons;
Mot1=zeros(6,1);
Mot1(3)=Param.tau1;
%Mot1(4:5)=F1;
%Mot2=zeros(6,1);
%Mot2(3)=Param.tau2;
%Mot1(4:5)=F2;
Fa=[Mot1;Act(1:Param.n*Param.na,:)];
res=M*qdd-Fa+Fe+Co*qd+Keps*q+Deps*qd; %2nd derivative of q
res=[res(3);res(7:end)];


 eps_q=Keps+Feq;
 eps_qd=Cqd+Deps;
% 
 eps_q=[eps_q(3,:);eps_q(7:end,:)];
 K_star=[eps_q(:,3) eps_q(:,7:end)];
% 
 eps_qd=[eps_qd(3,:);eps_qd(7:end,:)];
 C_star=[eps_qd(:,3) eps_qd(:,7:end)];
% 
 M=[M(3,:);M(7:end,:)];
 M=[M(:,3) M(:,7:end)];
 C=[C(3) C(7:9) C(12) C(16:end)];
% Psi_q=[Psi_q(:,3) Psi_q(:,7:9) Psi_q(:,12) Psi_q(:,16:end)];
% 
 M_inv=inv(M);
%  P=eye(2+Param.n*Param.na+3)-Psi_q'*inv(Psi_q*M_inv*Psi_q')*Psi_q*M_inv;
% 
%  M_star=M;
%  M_sinv=inv(M_star);
%  C_star=P*eps_qd;
% 
% K_star=P*eps_q;
H_star=-[1;zeros(Param.n*Param.na,1)]; 
% C=[zeros(2,2+Param.n*Param.na+1), eye(2) ,zeros(2,Param.n*Param.na+5)];
 A=[zeros(Param.n*Param.na+1), eye(Param.n*Param.na+1); -M_inv*K_star, -M_inv*C_star];
 B=[zeros(Param.n*Param.na+1,1);-M_inv*H_star];
sys = ss(A, B, C, 0);
G = tf(sys);
[V,eigen]= eig(A,"vector");
%[V1,eigen1]= eig(-M_sinv*K_star,"vector");
w=abs(eigen);

end