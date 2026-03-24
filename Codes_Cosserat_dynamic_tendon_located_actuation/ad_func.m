function a = ad_func(Y)
W=Y(1:3);
U=Y(4:6);
a=[Skew_symmetric(W) zeros(3) ; Skew_symmetric(U) Skew_symmetric(W)];
end