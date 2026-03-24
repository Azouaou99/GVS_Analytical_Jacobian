function H=Ftendons(X,i,q,Param) %Basis matrix
xi_a=Xi(X,q,Param);
xi=Param.B*xi_a + Param.B_bar*Param.xi_c; %Strain
K=xi(1:3); %angular strain
Gamma=xi(4:6);%Linear strain
H=zeros(6,4);
for m=1:4
Gamma_ci= Gamma + Skew_symmetric(K)*Param.Tendons_list(3*(i-1)+1:3*i,m);
H(:,m)=[Skew_symmetric(Param.Tendons_list(3*(i-1)+1:3*i,m))*Gamma_ci ; Gamma_ci]/norm(Gamma_ci);
end