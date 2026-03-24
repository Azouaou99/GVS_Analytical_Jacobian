function A = Ad_function(g)
R=g(1:3,1:3);
r=g(1:3,4);
A=[R zeros(3) ; Skew_symmetric(r)*R R];
end