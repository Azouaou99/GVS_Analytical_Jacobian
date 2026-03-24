%Main
%Parameters
clc
clear
format long
% Param.E = 10^6; % Young modulu
% Param.r = 0.01; % rod radius
% Param.g =0*[0 ; 0 ; 0 ; 0 ; 0; -9.81]; % gravity
% Param.L = 1; % Rod length
% Param.A =  pi*Param.r^2; % Cross section area
% Param.J1 = pi*Param.r^4/2; % Polar inertia moment
% Param.J2 = pi*Param.r^4/4; % Inertia moment
% Param.J3 = pi*Param.r^4/4; % Inertia moment
% %Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
% Param.H = diag([1.5708e-2 , 0.7854e-2 , 0.7854e-2 , 3.14e2 ,3.14e2 , 3.14e2]); %Hooke tensor

Param.E = 2.563*10^5; % Young modulu
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

Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);
%Param.mu=10^-3;
%Tendons parameters
%Param.Rb= 0.01/2; % Distance between a tendon and the backbone
Param.Rb= 0.02;
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)]; % coordinate of the tendon in a cross-section
Param.n=3;  %strain modes number
Param.na=2; %number of actuated strains
Param.DeltaX=Param.L/100;%length of segments
Param.DeltaX2=Param.L/100;
Param.n_seg=length(Param.DeltaX2:Param.DeltaX2:Param.L);
Param.n_segments=length(Param.DeltaX:Param.DeltaX:Param.L);%number of segments
Param.n_lambda=3; %dimension of Lambda
Y=0:Param.DeltaX2:Param.L;
Param.Tendons_list=zeros(3*(Param.n_seg+1),4);
for i=1:Param.n_seg+1
Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2, Param.Rb)]; 
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.DeltaX2)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.DeltaX2)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.DeltaX2)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.DeltaX2))]; % Convergent routing
%Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.DeltaX2)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.DeltaX2)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.DeltaX2)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.DeltaX2))]; % Convergent routing
%Parallel truncated
%if (i-1)*Param.DeltaX2<Param.L/2
%    Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)];
%end
end
%Param.Forces_Tendons = zeros(1,4); % tension in each tendon
   
   

  
Param.Forces_Tendons=[0;0;0;0];   
Param.Ftip=0*[0;-1/sqrt(2);1/sqrt(2);0;0;0];
%Param.Forces_Tendons(1)=1;


%Initial conditions
Param.g0=eye(4);
%Param.B =[eye(3);zeros(3)]; 
%Param.B_bar = [zeros(3);eye(3)];
% Param.B     = [0 0 0;1 0 0;0 0 0;0 1 0;0 0 0;0 0 1];
% Param.B_bar = [1 0 0;0 0 0;0 1 0;0 0 0;0 0 1;0 0 0];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0];

% Param.B     = [0;1;0;0;0;0];
% Param.B_bar = [1 0 0 0 0;0 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;1;0;0];

Param.B     = [0 0;1 0;0 0;0 1;0 0;0 0];
Param.B_bar = [1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0;0 0 1 0;0 0 0 1];
Param.xi_0=[0;0;0;1;0;0];
Param.xi_a0=Param.B'*Param.xi_0;
Param.xi_c=[0;0;0;0];

Param.Ha=Param.B'*Param.H*Param.B;%the matrix of the reduced Hooke coefficients
%Param.Da=Param.B'*Param.D*Param.B;%reduced damping matrix
Y=0:Param.DeltaX2:Param.L;
SK=zeros(Param.n*Param.na,Param.n*Param.na);
for k=1:Param.n_seg+1
    phival=Phi(Param.na,Param.n,Y(k),Param.L);
    fK=phival'*Param.Ha*phival;
    if k==1
      %  fK10=fK1;
        fK0=fK;
    elseif k==(Param.n_seg+1)
       % fK1n=fK1;
        fKn=fK;
    else
       % SK1=SK1+fK1;
        SK=SK+fK;
    end
end
Param.Keps= Param.DeltaX2.*((fK0+fKn)./2 + SK);
%Trapeze integration method
% init=[   6.429620303917535
%  -11.212951359843338
%    6.448899705153095
%    0.003766202295381
%   -0.001967763183432
%   -0.003079962039788];
init=zeros(Param.na*Param.n,1);
 for T=[0:0.01:5]
     Param.Forces_Tendons=[0;0;0;0]; 
     Param.Ftip=[0;T;0;0;0;0];
%Param.phip=zeros(4,1);
%Param.qp=ones(2*Param.na*Param.n,1);
%    %options = optimoptions("fsolve");
%      options=optimset('fsolve');
%      options.MaxFunEvals =1000000000000;
%      options.MaxIter = 50000000000;
%      options.TolFun=1e-5;
% %     % options.StepTolerance=10^-4;
% %     %options.exitflag=4;
% %     %options.OptimalityTolerance=10^(-2);
%     options.Jacobian='on';
%      options.Display = 'iter';
%     options.Algorithm='levenberg-marquardt';

    % options  =  optimoptions( 'fsolve', ...
    %                       'FiniteDifferenceStepSize', 1e-5, ...
    %                       'OptimalityTolerance', 1e-5, ...
    %                       'StepTolerance', 1e-5, ...
    %                       'FunctionTolerance', 1e-5, ...
    %                       'SpecifyObjectiveGradient', true, ...
    %                       'Display', 'off');
tic
%[q1,f,a,b,J1]=fsolve(@(var)Static_V2_grad(var,Param),init,options);
%J11=J1-Param.Keps
%options.Jacobian='off';
%[q1,f,a,b,J2]=fsolve(@(var)Static_V2(var,Param),init,options);
%J22=J2-Param.Keps
q1=Newton_r(init,1e-5,Param);
init=q1;
toc

r1=zeros(3,Param.n_X+1);
Q1=zeros(4,Param.n_X+1);
Q1(:,1)=[1;0;0;0];
for j=1:Param.n_X
  %   g=Geometric_model((j-1)*Param.dX,q,Param);
  %   r(:,j)=g(1:3,4);
        phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
        xia1=Param.xi_a0 + phi*q1;
        xi1=Param.B*xia1+Param.B_bar*Param.xi_c
        %xi(4:6)=[1;0;0];
        K1=xi1(1:3); %angular strain
        Gamma1=xi1(4:6);%Linear strain 
        R1=eye(3) + 2/(Q1(:,j)'*Q1(:,j)) * [-Q1(3,j)^2-Q1(4,j)^2, Q1(2,j)*Q1(3,j)-Q1(4,j)*Q1(1,j),Q1(2,j)*Q1(4,j) + Q1(3,j)*Q1(1,j) ; 
                                 Q1(2,j)*Q1(3,j)+Q1(4,j)*Q1(1,j), -Q1(2,j)^2-Q1(4,j)^2,Q1(3,j)*Q1(4,j) - Q1(2,j)*Q1(1,j) ; 
                                 Q1(2,j)*Q1(4,j)-Q1(3,j)*Q1(1,j), Q1(3,j)*Q1(4,j) + Q1(2,j)*Q1(1,j), -Q1(2,j)^2-Q1(3,j)^2];
        Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
                 K1(1), 0, K1(3), -K1(2);
                K1(2), -K1(3), 0, K1(1);
                K1(3), K1(2), -K1(1), 0 ] * Q1(:,j)/2;
        r_X1 = R1*Gamma1;
        Q1(:,j+1)=Q1(:,j) + Q_X1*Param.dX;
        r1(:,j+1)=r1(:,j) + r_X1*Param.dX;
        R1=eye(3) + 2/(Q1(:,j+1)'*Q1(:,j+1)) * [-Q1(3,j+1)^2-Q1(4,j+1)^2, Q1(2,j+1)*Q1(3,j+1)-Q1(4,j+1)*Q1(1,j+1),Q1(2,j+1)*Q1(4,j+1) + Q1(3,j+1)*Q1(1,j+1) ; 
                                 Q1(2,j+1)*Q1(3,j+1)+Q1(4,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(4,j+1)^2,Q1(3,j+1)*Q1(4,j+1) - Q1(2,j+1)*Q1(1,j+1) ; 
                                 Q1(2,j+1)*Q1(4,j+1)-Q1(3,j+1)*Q1(1,j+1), Q1(3,j+1)*Q1(4,j+1) + Q1(2,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(3,j+1)^2];
end  
    %hold on
    plot(r1(1,:),r1(3,:),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow
    %hold off
 end