function inv = Inverse_T(g)
R=g(1:3,1:3);
r=g(1:3,4);
inv=[R', -R'*r; 0 0 0 1];
end