function res = Dynamic(q1,q0,qd0,i,Param)

qd1=(q1-q0)/Param.dt;
qdd1=(qd1-qd0)/Param.dt;

SM1=zeros(Param.n*Param.na,Param.n*Param.na);
SC1=zeros(Param.n*Param.na,Param.n*Param.na);
SFe1=zeros(Param.n*Param.na,1);
SFa1=zeros(Param.n*Param.na,1);
g1=eye(4);
J1= zeros(6,Param.n*Param.na);
I1= zeros(6,Param.n*Param.na);
Jd1= zeros(6,Param.n*Param.na);
Id1= zeros(6,Param.n*Param.na);
Y=0:Param.dX:Param.L;
for k=1:Param.n_X+1
    %j=floor(Y(k)/Param.dX);
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    xiaval1=Param.xi_a0+phival*q1;
    xival1=Param.B*xiaval1+Param.B_bar*Param.xi_c;
    H1=Ftendons(Y(k),k,q1,Param);
    fFa1=-phival'*Param.B'*H1;
    fM1=J1'*Param.M*J1;
    fC1=J1'*(Param.M*Jd1 -ad_func(J1*qd1)' * Param.M*J1);
    fFe1=-J1'*Param.M*Ad_function(Inverse_T(g1))*Param.g;
    if k==1
        fM10=fM1;
        fC10=fC1;
        fFe10=fFe1;
        fFa10=fFa1;
    end
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
    if k~=Param.n_X+1
    P1=Ad_function(g1)*Param.B*phival;
    Pd1=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
    %g1=g1*expm(Hat(xival1)*(Param.dX));
    g1=g1*expm2(Hat(xival1)*(Param.dX),xival1*(Param.dX));
    N1=Ad_function(g1)*Param.B*phival;
    I1=I1+(P1+N1)*Param.dX/2;
    J1=Ad_function(Inverse_T(g1))*I1;
    Nd1=Ad_function(g1)*ad_func(Param.B*phival*qd1)*J1;
    Id1=Id1+(Pd1+Nd1)*Param.dX/2;
    Jd1=-Ad_function(Inverse_T(g1))*Id1;
    end
end
%model matrices projected on the basis function
Mq=Param.L/(Param.n_X)*((fM10+fM1n)/2 + SM1);
Cqqd=Param.L/(Param.n_X)*((fC10+fC1n)/2 + SC1);
Fe=Param.L/(Param.n_X)*((fFe10+fFe1n)/2 + SFe1);
Hq= Param.L/(Param.n_X)*((fFa10+fFa1n)/2 + SFa1);

res=Mq*qdd1-Hq*Param.Forces_Tendons+Fe+Cqqd*qd1+Param.Keps*q1+Param.Deps*qd1;%Residue
end