function res = Static(q,Param)
SFa1=zeros(Param.n*Param.na,1);
%Se1=zeros(Param.n*Param.na,1);
SFe1=zeros(Param.n*Param.na,1);
J1= zeros(6,Param.n*Param.na);
I1= zeros(6,Param.n*Param.na);
Y=0:Param.DeltaX2:Param.L;
g1=Param.g0;
for k=1:Param.n_seg+1
   % Jval1=Jacobian(Y(k),q1,Param);
   % Jval2=Jacobian(Y(k),q2,Param);
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    xiaval1=Param.xi_a0+phival*q;
    xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
%     K1=xival1(1:3); %angular strain
%     Gamma1=xival1(4:6);%Linear strain
%     R1=eye(3) + 2/(Q1'*Q1) * [-Q1(3)^2-Q1(4)^2, Q1(2)*Q1(3)-Q1(4)*Q1(1),Q1(2)*Q1(4) + Q1(3)*Q1(1) ;
%         Q1(2)*Q1(3)+Q1(4)*Q1(1), -Q1(2)^2-Q1(4)^2,Q1(3)*Q1(4) - Q1(2)*Q1(1) ;
%         Q1(2)*Q1(4)-Q1(3)*Q1(1), Q1(3)*Q1(4) + Q1(2)*Q1(1), -Q1(2)^2-Q1(3)^2];
%     Q1_X = [ 0, -K1(1), -K1(2), -K1(3);
%         K1(1), 0, K1(3), -K1(2);
%         K1(2), -K1(3), 0, K1(1);
%         K1(3), K1(2), -K1(1), 0 ] * Q1/2;
%     r1_X = R1*Gamma1;
    H1=Ftendons(Y(k),k,q,Param);
    fFa1=-phival'*Param.B'*H1;
    %fJ1=(Ad_function(g1))*Param.B*phival;
    fFe1=-J1'*Param.M*Ad_function(g1)^(-1)*Param.g;
    if k==1
       fFa10=fFa1;
        fFe10=fFe1;
    elseif k==(Param.n_seg+1)
        fFa1n=fFa1;
        fFe1n=fFe1;
    else
       SFa1=SFa1+fFa1;
       SFe1=SFe1+fFe1;
    end
    P1=Ad_function(g1)*Param.B*phival;
    %g1=g1*expm2(Hat(xival1)*(Param.DeltaX2),xival1*(Param.DeltaX2));
    g1=g1*expm(Hat(xival1)*(Param.DeltaX2));
    N1=Ad_function(g1)*Param.B*phival;
    I1=I1+(P1+N1)*Param.DeltaX2/2;
    J1=Ad_function(g1)^(-1)*I1;
%       Q1=Q1 + Q1_X*Param.DeltaX2;
%       r1=r1 + r1_X*Param.DeltaX2;
%       R1=eye(3) + 2/(Q1'*Q1) * [-Q1(3)^2-Q1(4)^2, Q1(2)*Q1(3)-Q1(4)*Q1(1),Q1(2)*Q1(4) + Q1(3)*Q1(1) ;
%       Q1(2)*Q1(3)+Q1(4)*Q1(1), -Q1(2)^2-Q1(4)^2,Q1(3)*Q1(4) - Q1(2)*Q1(1) ;
%       Q1(2)*Q1(4)-Q1(3)*Q1(1), Q1(3)*Q1(4) + Q1(2)*Q1(1), -Q1(2)^2-Q1(3)^2];
% 
%       Q2=Q2 + Q2_X*Param.DeltaX2;
%       r2=r2 + r2_X*Param.DeltaX2;
%       R2=eye(3) + 2/(Q2'*Q2) * [-Q2(3)^2-Q2(4)^2, Q2(2)*Q2(3)-Q2(4)*Q2(1),Q2(2)*Q2(4) + Q2(3)*Q2(1) ;
%         Q2(2)*Q2(3)+Q2(4)*Q2(1), -Q2(2)^2-Q2(4)^2,Q2(3)*Q2(4) - Q2(2)*Q2(1) ;
%         Q2(2)*Q2(4)-Q2(3)*Q2(1), Q2(3)*Q2(4) + Q2(2)*Q2(1), -Q2(2)^2-Q2(3)^2];
end
r1=g1(1:3,4);
R1=g1(1:3,1:3);
%[Jval1,Jval2]=Jacobian_2legs(Param.L,q1,q2,Param);
% Jval1=(Param.DeltaX2*((fJ10+fJ1n)./2 + SJ1)*(Ad_function(g1)^(-1)'))';
% Jval2=(Param.DeltaX2*((fJ20+fJ2n)./2 + SJ2)*(Ad_function(g2)^(-1)'))';
Fe=Param.DeltaX2*((fFe10+fFe1n)./2 + SFe1);
Hq= Param.DeltaX2.*((fFa10+fFa1n)./2 + SFa1);
Keps=Param.Keps;
res =Keps*q - Hq*Param.Forces_Tendons + Fe;
% Param.qp=q(1:2*Param.na*Param.n);
% Param.phip=phi;
end