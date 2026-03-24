function e=expm2(xi_h,xi) %Basis matrix
n_xi=norm(xi);
e=eye(4)+xi_h+1*norm(xi)^2*(1-cos(n_xi))*xi_h^2 + 1*n_xi^3*(n_xi-sin(n_xi))*xi_h^3;
end