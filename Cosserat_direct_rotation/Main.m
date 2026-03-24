%Main
%Parameters
clc
clear
format long
%  Param.E = 200e9; % Young modulu
%  Param.nu=0.3; % poisson's ratio
%  Param.G=Param.E/(2*(1+Param.nu)); % shear modulu
%  Param.r = 0.001; % rod radius
%  Param.rho = 8000; % Mass density
Param.E = 2.563*10^4; % Young modulu
Param.G=8.543*10^4; % shear modulu
Param.r = 0.01; % rod radius
Param.rho = 1.41*10^3; % Mass density

Param.g =0*[0 ; 0 ; 0 ; 9.81 ; 0; 0]; % gravity
Param.L = 0.2; % Rod length
Param.A = pi*Param.r^2; % Cross section area
Param.J1 = pi*Param.r^4/2; % Polar inertia moment
Param.J2 = pi*Param.r^4/4; % Inertia moment
Param.J3 = pi*Param.r^4/4; % Inertia moment
%Param.mu=1000;
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
%Param.D = Param.mu*diag([Param.J1 , 3*Param.J2 , 3*Param.J3 , 3*Param.A , Param.A , Param.A]);%Damping matrix
Param.D=10^(-4)*eye(6);%Damping matrix
Param.duree=4;%simulation time
Param.dt=0.01; %time step
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);
Param.n_t=length(Param.dt:Param.dt:Param.duree);
Param.n_lambda=6;

%Platform parameters
Param.rhop=0*1.41*10^3;
Param.a=0.02;
Param.b=0.08;
Param.c=0.02;
Param.Vp=Param.a*Param.b*Param.c;
Param.mp=Param.rhop*Param.Vp;
Param.Jp3=1/12*Param.mp*(Param.a^2+Param.b^2);
Param.Mp = diag([Param.Jp3,Param.mp,Param.mp]);%cross-sectional inertia
Param.F_p=[0;-Param.mp*9.81;0];
%Param.Fep=[0;20;0];

%Tendons parameters
Param.Rb= 0.025; % Distance between a tendon and the backbone
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)];% coordinate of the tendon in a cross-section
%Param.Tendons_list = [Param.Tendon_coordinate(0),Param.Tendon_coordinate(pi/2),Param.Tendon_coordinate(pi),Param.Tendon_coordinate(3*pi/2)]; % 4 tendons separated with an angle of pi/2
Param.n=3;  %strain modes number
Param.na=2; %number of actuated strains

Y=0:Param.dX:Param.L;
Param.Tendons_list=zeros(3*(Param.n_X+1),2);
for i=1:Param.n_X+1
Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb)]; % Parallel case
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.dX))]; % Convergent routing
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.dX))]; % Convergent routing
%Parallel truncated
%if (i-1)*Param.dX<Param.L/2
%    Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)];
%end
end
%Param.Forces_Tendons = zeros(4,1); % tension in each tendon
   

%  Param.Forces_Tendons=[0
%     16.426877086020717
%     18.192266697163880
% 0]; 

  Param.Forces_Tendons=[0
    1
    0
0]; 
%Param.Ftip=[0;0;0;0;0;2];

%Initial conditions
Param.c1=0*pi/500;
Param.c2=0*-pi/500;
Param.theta01=0*-pi/8;
Param.theta02=0*pi/8;
R01=[cos(Param.theta01) -sin(Param.theta01) 0; sin(Param.theta01) cos(Param.theta01) 0; 0 0 1];
R02=[cos(Param.theta02) -sin(Param.theta02) 0; sin(Param.theta02) cos(Param.theta02) 0; 0 0 1];
Param.g1=eye(4);
Param.g1(1:3,1:3)=R01;
Param.g2=eye(4);
Param.y2=0.034;
Param.g2(1:3,1:3)=R02;
Param.g2(2,4)=Param.y2;
Param.CM1=[Param.a/2;Param.y2/2];
Param.CM2=[Param.a/2;-Param.y2/2];

% Param.B =[0 0;1 0;0 0;0 0;0 0;0 1];
% Param.B_bar = [1 0 0 0;0 0 0 0; 0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0;1];

Param.B     = [0 0;0 0;1 0;0 1;0 0;0 0];%Selection matrix
Param.B_bar = [1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0;0 0 1 0;0 0 0 1];
Param.xi_0=[0;0;0;1;0;0];%reference configuration
Param.xi_a0=Param.B'*Param.xi_0;
Param.xi_c=[0;0;0;0];

% Param.B     = [0;0;1;0;0;0];%Selection matrix
% Param.B_bar = [1 0 0 0 0;0 1 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
% Param.xi_0=[0;0;0;1;0;0];%reference configuration
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;1;0;0];


% Param.B     = [0 0 0 ;1 0 0;0 0 0;0 1 0;0 0 0;0 0 1];
% Param.B_bar = [1 0 0;0 0 0;0 1 0;0 0 0;0 0 1;0 0 0];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0];

Param.Ha=Param.B'*Param.H*Param.B;%the matrix of the reduced Hooke coefficients
Param.Da=Param.B'*Param.D*Param.B;%reduced damping matrix

  
  
   


%Initial conditions
% Param.q10=[3.67479702448598e-05
% 0.000296870329238203
% 0.000133714523684237
% 0.0431893381573064
% -0.00669064246093489
% 0.000269545619105613];
% 
% Param.q20=[3.67466172407047e-05
% 0.000296733783306033
% 0.000133749039599706
% 0.0431899206335046
% -0.00669065412636491
% -0.00163709037943777];
% 
% Param.qp0=[6.75623552187804e-06
% 0.228651393425217
% 0.0249983768856270];
  
   
 Param.q10=zeros(Param.na*Param.n,1);
  % Param.q10=[ 4.366271598191997
  % -0.000000000000501
  %  5.739081835533309
  % -0.184832887744930
  %  0.000000000000005
  % -0.000998259329265
  %  ];
% 
Param.q20=zeros(Param.na*Param.n,1);
 %   Param.q20=[   4.366271597999246
 %  -0.000000000007372
 % -17.607626210206082
 %  -0.243326785302555
 %  -0.000000000000010
 %   0.001891688868647];
% 
  Param.qp0=[0;0.2+Param.a/2;Param.y2/2];
     % Param.qp0=[ 0.873254319625532
     % 0.14
     % 0.085];

 Param.q0=[Param.q10;Param.q20;Param.qp0];
 %[3.926988884408633
%   -0.000187366956875
%   -0.005875132022834
%    0.000000000087274
%   -0.000002887195828
%    0.000000144241633
%   -3.926990178091315
%    0.000187287566010
%    0.005875575626324
%    0.000000000373692
%   -0.000002887190550
%   -0.000004808138983
%   -0.000000129284439
%    0.204896040335919
%    0.024999984497435];
   
  
  


Param.qdp=zeros(3,1);
Param.qd0=zeros(Param.n*Param.na,1);
%initialization of the coefficients
q1=zeros(Param.na*Param.n,Param.n_t+1);
q2=zeros(Param.na*Param.n,Param.n_t+1);
qp=zeros(3,Param.n_t+1);
q=zeros(2*Param.na*Param.n+3,Param.n_t+1);

q1(:,1)=Param.q10;
q2(:,1)=Param.q20;
qp(:,1)=Param.qp0;

q(:,1)=Param.q0;

qd1=zeros(Param.na*Param.n,Param.n_t+1);
qd2=zeros(Param.na*Param.n,Param.n_t+1);
qdp=zeros(3,Param.n_t+1);
qd1(:,1)=Param.qd0;
qd2(:,1)=Param.qd0;
qdp(:,1)=Param.qdp;
qd=zeros(2*Param.na*Param.n+3,Param.n_t+1);


  lambda=zeros(Param.n_lambda,Param.n_t+1);
  % lambda(:,1)= [ 0.216500311897884
  %  1.425171903665252
  %  0.665080221710548
  % -0.237877890452878
  % -1.425171903665252
  % -0.665080221710548];
qdd1=zeros(Param.na*Param.n,Param.n_t+1);
qdd2=zeros(Param.na*Param.n,Param.n_t+1);
qddp=zeros(3,Param.n_t+1);


%r=zeros(3,Param.n_X+1,Param.n_t+1);
%Q=zeros(4,Param.n_X+1);
%Q(:,1)=[1;0;0;0];

r1=zeros(3,Param.n_X+1,Param.n_t+1);
r2=zeros(3,Param.n_X+1,Param.n_t+1);
r2(2,1,:)=Param.y2;
% for j=1:Param.n_X
%     phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
%     xia=Param.xi_a0 + phi*q(:,1);
%     xi=Param.B*xia+Param.B_bar*Param.xi_c;
%     K=xi(1:3); %angular strain
%     Gamma=xi(4:6);%Linear strain
%     R=eye(3) + 2/(Q(:,j)'*Q(:,j)) * [-Q(3,j)^2-Q(4,j)^2, Q(2,j)*Q(3,j)-Q(4,j)*Q(1,j),Q(2,j)*Q(4,j) + Q(3,j)*Q(1,j) ;
%         Q(2,j)*Q(3,j)+Q(4,j)*Q(1,j), -Q(2,j)^2-Q(4,j)^2,Q(3,j)*Q(4,j) - Q(2,j)*Q(1,j) ;
%         Q(2,j)*Q(4,j)-Q(3,j)*Q(1,j), Q(3,j)*Q(4,j) + Q(2,j)*Q(1,j), -Q(2,j)^2-Q(3,j)^2];
%     Q_X = [ 0, -K(1), -K(2), -K(3);
%         K(1), 0, K(3), -K(2);
%         K(2), -K(3), 0, K(1);
%         K(3), K(2), -K(1), 0 ] * Q(:,j)/2;
%     r_X = R*Gamma;
%     Q(:,j+1)=Q(:,j) + Q_X*Param.dX;
%     r(:,j+1,1)=r(:,j,1) + r_X*Param.dX;
% end
Q1=zeros(4,Param.n_X+1);
Q2=zeros(4,Param.n_X+1);
Q1(:,1)=[cos(Param.theta01/2);0;0;sin(Param.theta01/2)];
Q2(:,1)=[cos(Param.theta02/2);0;0;sin(Param.theta02/2)];
for j=1:Param.n_X
    phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
    xia1=Param.xi_a0 + phi*q1(:,1);
    xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
    xia2=Param.xi_a0 + phi*q2(:,1);
    xi2=Param.B*xia2+Param.B_bar*Param.xi_c;
    %xi(4:6)=[1;0;0];
    K1=xi1(1:3); %angular strain
    Gamma1=xi1(4:6);%Linear strain
    K2=xi2(1:3); %angular strain
    Gamma2=xi2(4:6);%Linear strain

    R1=eye(3) + 2/(Q1(:,j)'*Q1(:,j)) * [-Q1(3,j)^2-Q1(4,j)^2, Q1(2,j)*Q1(3,j)-Q1(4,j)*Q1(1,j),Q1(2,j)*Q1(4,j) + Q1(3,j)*Q1(1,j) ;
        Q1(2,j)*Q1(3,j)+Q1(4,j)*Q1(1,j), -Q1(2,j)^2-Q1(4,j)^2,Q1(3,j)*Q1(4,j) - Q1(2,j)*Q1(1,j) ;
        Q1(2,j)*Q1(4,j)-Q1(3,j)*Q1(1,j), Q1(3,j)*Q1(4,j) + Q1(2,j)*Q1(1,j), -Q1(2,j)^2-Q1(3,j)^2];
    Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
        K1(1), 0, K1(3), -K1(2);
        K1(2), -K1(3), 0, K1(1);
        K1(3), K1(2), -K1(1), 0 ] * Q1(:,j)/2;
    r_X1 = R1*Gamma1;
    Q1(:,j+1)=Q1(:,j) + Q_X1*Param.dX;
    r1(:,j+1,1)=r1(:,j,1) + r_X1*Param.dX;
    %  R1=eye(3) + 2/(Q1(:,j+1)'*Q1(:,j+1)) * [-Q1(3,j+1)^2-Q1(4,j+1)^2, Q1(2,j+1)*Q1(3,j+1)-Q1(4,j+1)*Q1(1,j+1),Q1(2,j+1)*Q1(4,j+1) + Q1(3,j+1)*Q1(1,j+1) ;
    %                           Q1(2,j+1)*Q1(3,j+1)+Q1(4,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(4,j+1)^2,Q1(3,j+1)*Q1(4,j+1) - Q1(2,j+1)*Q1(1,j+1) ;
    %                           Q1(2,j+1)*Q1(4,j+1)-Q1(3,j+1)*Q1(1,j+1), Q1(3,j+1)*Q1(4,j+1) + Q1(2,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(3,j+1)^2];


    R2=eye(3) + 2/(Q2(:,j)'*Q2(:,j)) * [-Q2(3,j)^2-Q2(4,j)^2, Q2(2,j)*Q2(3,j)-Q2(4,j)*Q2(1,j),Q2(2,j)*Q2(4,j) + Q2(3,j)*Q2(1,j) ;
        Q2(2,j)*Q2(3,j)+Q2(4,j)*Q2(1,j), -Q2(2,j)^2-Q2(4,j)^2,Q2(3,j)*Q2(4,j) - Q2(2,j)*Q2(1,j) ;
        Q2(2,j)*Q2(4,j)-Q2(3,j)*Q2(1,j), Q2(3,j)*Q2(4,j) + Q2(2,j)*Q2(1,j), -Q2(2,j)^2-Q2(3,j)^2];
    Q_X2 = [ 0, -K2(1), -K2(2), -K2(3);
        K2(1), 0, K2(3), -K2(2);
        K2(2), -K2(3), 0, K2(1);
        K2(3), K2(2), -K2(1), 0 ] * Q2(:,j)/2;
    r_X2 = R2*Gamma2;
    Q2(:,j+1)=Q2(:,j) + Q_X2*Param.dX;
    r2(:,j+1,1)=r2(:,j,1) + r_X2*Param.dX;
    %R2=eye(3) + 2/(Q2(:,j+1)'*Q2(:,j+1)) * [-Q2(3,j+1)^2-Q2(4,j+1)^2, Q2(2,j+1)*Q2(3,j+1)-Q2(4,j+1)*Q2(1,j+1),Q2(2,j+1)*Q2(4,j+1) + Q2(3,j+1)*Q2(1,j+1) ;
    %                         Q2(2,j+1)*Q2(3,j+1)+Q2(4,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(4,j+1)^2,Q2(3,j+1)*Q2(4,j+1) - Q2(2,j+1)*Q2(1,j+1) ;
    %                         Q2(2,j+1)*Q2(4,j+1)-Q2(3,j+1)*Q2(1,j+1), Q2(3,j+1)*Q2(4,j+1) + Q2(2,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(3,j+1)^2];
end

Y=0:Param.dX:Param.L;
SK=zeros(Param.n*Param.na,Param.n*Param.na);
SD=zeros(Param.n*Param.na,Param.n*Param.na);
for k=1:Param.n_X+1
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    fK=phival'*Param.Ha*phival;
    fD=phival'*Param.Da*phival;
    if k==1
        fK0=fK;
        fD0=fD;
    end
    if k==Param.n_X+1
        fKn=fK;
        fDn=fD;
    end
    if and(k~=Param.n_X+1,k~=1)
        SK=SK+fK;
        SD=SD+fD;
    end
end
Param.Keps= Param.L/(Param.n_X)*((fK0+fKn)/2 + SK);%Stiffness matrix on the modal space
Param.Deps= Param.L/(Param.n_X)*((fD0+fDn)/2 + SD);%Damping matrix on the modal space

Param.Y=0:Param.dX:Param.L;
l=1;
Param.F_p=[0;0;0];
T=1;
   for i=1:Param.n_t %Time updating
       %T=T+0.05;
         Param.Forces_Tendons=[0
    T
    0
0]; 
   % Param.Fp=Param.Fp+[0;1;0];    
         % if i==20
         %     Param.F_p=[0;1;0];
         % end
%       Param.Forces_Tendons=[0
%  50
%  50
%  0]; 
%  else
%        Param.Forces_Tendons=[0
% 51
% 50
% 0]; 
%  end
%       if i>Param.n_t/2 & i<Param.n_t/2+10
%  Param.Ftendons=[1;1;0;1];
%       else 
%        Param.Ftendons=[1;0;0;1];   
%      end
tic
    disp (i*Param.dt)
    init=[q1(:,i);q2(:,i);qp(:,i);lambda(:,i)];
   %options = optimoptions("fsolve");
%      options=optimset('fsolve');
%      options.MaxFunEvals =1000000000000;
%      options.MaxIter = 50000000000;
%      options.TolFun=1e-5;
% %     % options.StepTolerance=10^-4;
% %     %options.exitflag=4;
% %     %options.OptimalityTolerance=10^(-2);
%     options.Jacobian='on';
%      options.Display = 'iter';
    % %options.Algorithm='levenberg-marquardt';

    options  =  optimoptions( 'fsolve', ...
                          'FiniteDifferenceStepSize', 1e-7, ...
                          'OptimalityTolerance', 1e-7, ...
                          'StepTolerance', 1e-7, ...
                          'FunctionTolerance', 1e-7, ...
                          'SpecifyObjectiveGradient', true, ...
                          'Display', 'off');
    x=fsolve(@(var)Dynamic_grad(var,q(:,i),qd(:,i),i,Param),init,options);%2nd derivative of q
    toc
    Param.Fp=[0;0;0];
    %Param.Forces_Tendons(1,1)=Param.Forces_Tendons(1,1)+0.01;
   % Param.Forces_Tendons(1,2)=Param.Forces_Tendons(1,2)+0.01;
    %Param.Forces_Tendons(2,3)=Param.Forces_Tendons(2,3)+0.01;
   % Param.Forces_Tendons(2,4)=Param.Forces_Tendons(2,4)+0.01;
%     qd1(:,i+1)=x(1:Param.n*Param.na);
%     qd2(:,i+1)=x(Param.n*Param.na+1:2*Param.n*Param.na);
%     qd(:,i+1)=x(1:2*Param.n*Param.na);
    q1(:,i+1)=x(1:Param.n*Param.na);
    q2(:,i+1)=x(Param.n*Param.na+1:2*Param.n*Param.na);
    qp(:,i+1)=x(2*Param.n*Param.na+1:2*Param.n*Param.na+3);
    q(:,i+1)=x(1:2*Param.n*Param.na+3);
    lambda(:,i+1)=x(2*Param.n*Param.na+4 : 2*Param.n*Param.na+3+Param.n_lambda);

    qd1(:,i+1)=(q1(:,i+1)-q1(:,i))/Param.dt;%x(2*Param.n*Param.na+1:3*Param.n*Param.na);%time intgeration
    qd2(:,i+1)=(q2(:,i+1)-q2(:,i))/Param.dt;%x(3*Param.n*Param.na+1:4*Param.n*Param.na);
    qdp(:,i+1)=(qp(:,i+1)-qp(:,i))/Param.dt;%x(3*Param.n*Param.na+1:4*Param.n*Param.na);
    qd(:,i+1)=(q(:,i+1)-q(:,i))/Param.dt;

    Q1(:,1)=[cos(Param.theta01/2);0;0;sin(Param.theta01/2)];
    Q2(:,1)=[cos(Param.theta02/2);0;0;sin(Param.theta02/2)];
    for j=1:Param.n_X
        phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
        xia1=Param.xi_a0 + phi*q1(:,i+1);
        xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
        xia2=Param.xi_a0 + phi*q2(:,i+1);
        xi2=Param.B*xia2+Param.B_bar*Param.xi_c;
        %xi(4:6)=[1;0;0];
        K1=xi1(1:3); %angular strain
        Gamma1=xi1(4:6);%Linear strain
        K2=xi2(1:3); %angular strain
        Gamma2=xi2(4:6);%Linear strain

        R1=eye(3) + 2/(Q1(:,j)'*Q1(:,j)) * [-Q1(3,j)^2-Q1(4,j)^2, Q1(2,j)*Q1(3,j)-Q1(4,j)*Q1(1,j),Q1(2,j)*Q1(4,j) + Q1(3,j)*Q1(1,j) ;
            Q1(2,j)*Q1(3,j)+Q1(4,j)*Q1(1,j), -Q1(2,j)^2-Q1(4,j)^2,Q1(3,j)*Q1(4,j) - Q1(2,j)*Q1(1,j) ;
            Q1(2,j)*Q1(4,j)-Q1(3,j)*Q1(1,j), Q1(3,j)*Q1(4,j) + Q1(2,j)*Q1(1,j), -Q1(2,j)^2-Q1(3,j)^2];
        Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
            K1(1), 0, K1(3), -K1(2);
            K1(2), -K1(3), 0, K1(1);
            K1(3), K1(2), -K1(1), 0 ] * Q1(:,j)/2;
        r_X1 = R1*Gamma1;
        Q1(:,j+1)=Q1(:,j) + Q_X1*Param.dX;
        r1(:,j+1,i+1)=r1(:,j,i+1) + r_X1*Param.dX;
        %  R1=eye(3) + 2/(Q1(:,j+1)'*Q1(:,j+1)) * [-Q1(3,j+1)^2-Q1(4,j+1)^2, Q1(2,j+1)*Q1(3,j+1)-Q1(4,j+1)*Q1(1,j+1),Q1(2,j+1)*Q1(4,j+1) + Q1(3,j+1)*Q1(1,j+1) ;
        %                           Q1(2,j+1)*Q1(3,j+1)+Q1(4,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(4,j+1)^2,Q1(3,j+1)*Q1(4,j+1) - Q1(2,j+1)*Q1(1,j+1) ;
        %                           Q1(2,j+1)*Q1(4,j+1)-Q1(3,j+1)*Q1(1,j+1), Q1(3,j+1)*Q1(4,j+1) + Q1(2,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(3,j+1)^2];


        R2=eye(3) + 2/(Q2(:,j)'*Q2(:,j)) * [-Q2(3,j)^2-Q2(4,j)^2, Q2(2,j)*Q2(3,j)-Q2(4,j)*Q2(1,j),Q2(2,j)*Q2(4,j) + Q2(3,j)*Q2(1,j) ;
            Q2(2,j)*Q2(3,j)+Q2(4,j)*Q2(1,j), -Q2(2,j)^2-Q2(4,j)^2,Q2(3,j)*Q2(4,j) - Q2(2,j)*Q2(1,j) ;
            Q2(2,j)*Q2(4,j)-Q2(3,j)*Q2(1,j), Q2(3,j)*Q2(4,j) + Q2(2,j)*Q2(1,j), -Q2(2,j)^2-Q2(3,j)^2];
        Q_X2 = [ 0, -K2(1), -K2(2), -K2(3);
            K2(1), 0, K2(3), -K2(2);
            K2(2), -K2(3), 0, K2(1);
            K2(3), K2(2), -K2(1), 0 ] * Q2(:,j)/2;
        r_X2 = R2*Gamma2;
        Q2(:,j+1)=Q2(:,j) + Q_X2*Param.dX;
        r2(:,j+1,i+1)=r2(:,j,i+1) + r_X2*Param.dX;
        %R2=eye(3) + 2/(Q2(:,j+1)'*Q2(:,j+1)) * [-Q2(3,j+1)^2-Q2(4,j+1)^2, Q2(2,j+1)*Q2(3,j+1)-Q2(4,j+1)*Q2(1,j+1),Q2(2,j+1)*Q2(4,j+1) + Q2(3,j+1)*Q2(1,j+1) ;
        %                         Q2(2,j+1)*Q2(3,j+1)+Q2(4,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(4,j+1)^2,Q2(3,j+1)*Q2(4,j+1) - Q2(2,j+1)*Q2(1,j+1) ;
        %                         Q2(2,j+1)*Q2(4,j+1)-Q2(3,j+1)*Q2(1,j+1), Q2(3,j+1)*Q2(4,j+1) + Q2(2,j+1)*Q2(1,j+1), -Q2(2,j+1)^2-Q2(3,j+1)^2];
    end


% l_rec=0.02;
% for j=1:Param.n_t+1
%    % j=50
%     plot3(r1(1,:,j),r1(2,:,j),r1(3,:,j),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%     xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
%     hold on
%     plot3(r2(1,:,j),r2(2,:,j),r2(3,:,j),'b','LineWidth', 4); title('Cosserat rod');drawnow
%     %xrec=[r1(1,end,j);r1(1,end,j)+R1(1,1)*l_rec ; r2(1,end)+R2(1,1)*l_rec ; r2(1,end) ; r1(1,end)];
%     % yrec=[r1(2,end);r1(2,end)+R1(2,1)*l_rec;r2(2,end)+R2(2,1)*l_rec;r2(2,end);r1(2,end)];
%     % plot(xrec,yrec,'r','LineWidth', 4)
%     hold off
% end
l_rec=0.02;
%for j=1:Param.n_t+1
   % j=50
    plot(r1(1,:,i),r1(2,:,i),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    hold on
    plot(r2(1,:,i),r2(2,:,i),'b','LineWidth', 4); title('Cosserat rod');drawnow
    plot(qp(2,i),qp(3,i),'g*','LineWidth', 4); title('Cosserat rod');drawnow
     xrec=[qp(2,i)+(Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i))); qp(2,i)+(Param.a/2*cos(qp(1,i))+Param.b/2*sin(qp(1,i))) ;qp(2,i)-(Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i)));qp(2,i)+(-Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i)));qp(2,i)+(Param.a/2*cos(qp(1,i))-Param.b/2*sin(qp(1,i)))];
     yrec=[qp(3,i)+(+Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i))) ;qp(3,i)+(Param.a/2*sin(qp(1,i))-Param.b/2*cos(qp(1,i))) ;qp(3,i)-(Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i)));qp(3,i)+(-Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i)));qp(3,i)+(Param.a/2*sin(qp(1,i))+Param.b/2*cos(qp(1,i)))];
     plot(xrec,yrec,'r','LineWidth', 4);drawnow
    %xrec=[r1(1,end,j);r1(1,end,j)+R1(1,1)*l_rec ; r2(1,end)+R2(1,1)*l_rec ; r2(1,end) ; r1(1,end)];
    % yrec=[r1(2,end);r1(2,end)+R1(2,1)*l_rec;r2(2,end)+R2(2,1)*l_rec;r2(2,end);r1(2,end)];
    % plot(xrec,yrec,'r','LineWidth', 4)
    hold off
end
 j=1;
 plot(r1(1,:,j),r1(2,:,j),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    hold on
    plot(r2(1,:,j),r2(2,:,j),'b','LineWidth', 4); title('Cosserat rod');drawnow
    plot(qp(2,j),qp(3,j),'g*','LineWidth', 4); title('Cosserat rod');;drawnow
     xrec=[qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j))); qp(2,j)+(Param.a/2*cos(qp(1,j))+Param.b/2*sin(qp(1,j))) ;qp(2,j)-(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(-Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)))];
     yrec=[qp(3,j)+(+Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j))) ;qp(3,j)+(Param.a/2*sin(qp(1,j))-Param.b/2*cos(qp(1,j))) ;qp(3,j)-(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(-Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)))];
     plot(xrec,yrec,'r','LineWidth', 4);drawnow
im={}; 
[im{j},map]=frame2im(getframe);
for j=2:Param.n_t+1
    figure(1)
    clf
    disp(j);
   plot(r1(1,:,j),r1(2,:,j),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    hold on
    plot(r2(1,:,j),r2(2,:,j),'b','LineWidth', 4); title('Cosserat rod');drawnow
    plot(qp(2,j),qp(3,j),'g*','LineWidth', 4); title('Cosserat rod');;drawnow
     xrec=[qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j))); qp(2,j)+(Param.a/2*cos(qp(1,j))+Param.b/2*sin(qp(1,j))) ;qp(2,j)-(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(-Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)));qp(2,j)+(Param.a/2*cos(qp(1,j))-Param.b/2*sin(qp(1,j)))];
     yrec=[qp(3,j)+(+Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j))) ;qp(3,j)+(Param.a/2*sin(qp(1,j))-Param.b/2*cos(qp(1,j))) ;qp(3,j)-(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(-Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)));qp(3,j)+(Param.a/2*sin(qp(1,j))+Param.b/2*cos(qp(1,j)))];
     plot(xrec,yrec,'r','LineWidth', 4);drawnow
    %xrec=[r1(1,end,j);r1(1,end,j)+R1(1,1)*l_rec ; r2(1,end)+R2(1,1)*l_rec ; r2(1,end) ; r1(1,end)];
    % yrec=[r1(2,end);r1(2,end)+R1(2,1)*l_rec;r2(2,end)+R2(2,1)*l_rec;r2(2,end);r1(2,end)];
%     % plot(xrec,yrec,'r','LineWidth', 4)
    [im{j},map]=frame2im(getframe);
hold off;
 pause(0.001)
end
[temp,map]=rgb2ind(im{1},Param.n_t+1);
for j=1:Param.n_t+1
    gifim(:,:,1,j)=rgb2ind(im{j},map);
end
imwrite(gifim,map,'Stability.gif');


l_ana=zeros(Param.n_t+1,1);
m_ana=zeros(Param.n_t+1,1);
n_ana=zeros(Param.n_t+1,1);
t=(0:Param.dt:Param.duree)';
for j=1:Param.n_t+1
    l_ana(j,1)=(r1(1,end,j)+r2(1,end,j))/2;
    m_ana(j,1)=(r1(2,end,j)+r2(2,end,j))/2;
    n_ana(j,1)=(r1(3,end,j)+r2(3,end,j))/2;
end
figure; hold on
grid
plot(t',l_ana,'r')
plot(t',m_ana,'b')
plot(t',n_ana,'g')
ylabel('$x$, $y$, $z$ [m] ','interpreter', 'latex');
xlabel('$t$ [s]','interpreter', 'latex');
legend({'$x_{eff}$', '$y_{eff}$', '$z_{eff}$'}, 'interpreter', 'latex')
title('End effector positions')

% figure; hold on
% grid
% plot(t',qp(2,:),'r')
% plot(t',qp(3,:),'b')
% %plot(t',n_ana,'g')
% ylabel('$x$, $y$ [m] ','interpreter', 'latex');
% xlabel('$t$ [s]','interpreter', 'latex');
% legend({'$x_{p}$', '$y_{p}$'}, 'interpreter', 'latex')
% title('End effector positions')


