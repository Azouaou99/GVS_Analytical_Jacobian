function g = Geometric_model(X,q,Param)%geometric model integration
j=floor(X/Param.DeltaX);
xia_j=Param.xi_a0 + Phi(Param.na,Param.n,j*Param.DeltaX,Param.L)*q;
xi_j=Param.B*xia_j+Param.B_bar*Param.xi_c;
g=Param.g0;
if j~=0
    for i=0:j-1
        xia_i=Param.xi_a0 + Phi(Param.na,Param.n,i*Param.DeltaX,Param.L)*q;
        xi_i=Param.B*xia_i+Param.B_bar*Param.xi_c;
        g=g*expm(Hat(xi_i)*Param.DeltaX);
    end
end
g=g*expm(Hat(xi_j)*(X-j*Param.DeltaX));
end
