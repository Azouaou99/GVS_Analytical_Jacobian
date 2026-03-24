function q1=Newton_r(q0,qp,qdp,err,Param)
i=0;
f=10;
while and(max(f)>err,i<15)
   [f,grad]= Dynamic_grad_NR(q0,qp,qdp,Param);
    q1=q0-grad\f;
    i=i+1;
    q0=q1;
end
end