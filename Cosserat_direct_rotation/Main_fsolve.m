%Main
%Parameters
clc
clear
tic
format long
%% Mechanical characteristics
Param.E = 2.563*10^6; % Young modulu
Param.G=8.543*10^4; % shear modulu
Param.r = 0.01; % rod radius
Param.rho = 1.41*10^3; % Mass density
Param.g =0*[0 ; 0 ; 0 ; -9.81 ; 0; 0]; % gravity
Param.L = 0.2; % Rod length
Param.A = pi*Param.r^2; % Cross section area
Param.J1 = pi*Param.r^4/2; % Polar inertia moment
Param.J2 = pi*Param.r^4/4; % Inertia moment
Param.J3 = pi*Param.r^4/4; % Inertia moment
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
Param.D=0*10^(-6)*eye(6);%Damping matrix
Param.duree=5;%simulation time
Param.dt=0.001; %time step
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);
Param.n_t=length(Param.dt:Param.dt:Param.duree);
Param.n=3;  %strain modes number
Param.na=1; %number of actuated strains

%% Tendons parameters
Param.Rb= 0.02; % Distance between a tendon and the backbone
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)];% coordinate of the tendon in a cross-section
Param.Tendons_list=zeros(3*(Param.n_X+1),2);
for i=1:Param.n_X+1
    % Parallel case
    Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb)]; 
    
    % Convergent routing
    %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.dX))]; 
    
    %Cross routing
    %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.dX))]; 

    %Parallel truncated
    % if (i-1)*Param.dX<Param.L/2
    %     Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)];
    % end
end
Param.Forces_Tendons=[0;0]; % tension in each tendon
Param.tau1=0;%-0*sin(40*i*Param.dt);
%% DOF Reduction

%bending in y and extension in x
% Param.B     = [0 0;1 0;0 0;0 1;0 0;0 0];%Selection matrix
% Param.B_bar = [1 0 0 0;0 0 0 0;0 1 0 0;0 0 0 0;0 0 1 0;0 0 0 1];
% Param.xi_0=[0;0;0;1;0;0];%reference configuration
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0;0];

%bending in y, extension in x and shear in z
% Param.B     = [0 0 0 ;1 0 0;0 0 0;0 1 0;0 0 0;0 0 1];
% Param.B_bar = [1 0 0;0 0 0;0 1 0;0 0 0;0 0 1;0 0 0];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=[0;0;0];

Param.B     = [0;0;1;0;0;0];%Selection matrix
Param.B_bar = [1 0 0 0 0;0 1 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
Param.xi_0=[0;0;0;1;0;0];%reference configuration
Param.xi_a0=Param.B'*Param.xi_0;
Param.xi_c=[0;0;1;0;0];

% Param.B     = [0 0 0 ;0 1 0;0 0 1;1 0 0;0 0 0;0 0 0];
% Param.B_bar = [1 0 0;0 0 0;0 0 0;0 0 0;0 1 0;0 0 1];
% Param.xi_0=[0;0;0;1;0;0];
% Param.xi_a0=Param.B'*Param.xi_0;
% Param.xi_c=Param.B_bar'*Param.xi_0;

Param.Ha=Param.B'*Param.H*Param.B;%the matrix of the reduced Hooke coefficients
Param.Da=Param.B'*Param.D*Param.B;%reduced damping matrix

%% Initialization
q1=zeros(Param.na*Param.n,Param.n_t+1);%Generalized coefficients
qa=zeros(1,Param.n_t+1);
qa_d6=zeros(6,Param.n_t+1);
qa_d6(3,:)=qa(:);
qd1=zeros(Param.na*Param.n,Param.n_t+1);%Generalized velocities
eta0=zeros(6,Param.n_t+1);
r1=zeros(3,Param.n_X+1,Param.n_t+1);%Positions of the robot sections
Q1=zeros(4,Param.n_X+1); %Quaternions
Q1(:,1)=[cos(qa(1)/2);0;0;sin(qa(1)/2)];
R01=[cos(qa(1)) -sin(qa(1)) 0; sin(qa(1)) cos(qa(1)) 0; 0 0 1];
Param.g0=eye(4);
Param.g0(1:3,1:3)=R01;
%% Beam reconstruction at t=0
for j=1:Param.n_X
    phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
    xia1=Param.xi_a0 + phi*q1(:,1);
    xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
    %xi(4:6)=[1;0;0];
    K1=xi1(1:3); %angular strain
    Gamma1=xi1(4:6);%Linear strain

    R1=eye(3) + 2/(Q1(:,j)'*Q1(:,j)) * [-Q1(3,j)^2-Q1(4,j)^2, Q1(2,j)*Q1(3,j)-Q1(4,j)*Q1(1,j),Q1(2,j)*Q1(4,j) + Q1(3,j)*Q1(1,j) ; %Relation between quaternions and rotation matrix
        Q1(2,j)*Q1(3,j)+Q1(4,j)*Q1(1,j), -Q1(2,j)^2-Q1(4,j)^2,Q1(3,j)*Q1(4,j) - Q1(2,j)*Q1(1,j) ;
        Q1(2,j)*Q1(4,j)-Q1(3,j)*Q1(1,j), Q1(3,j)*Q1(4,j) + Q1(2,j)*Q1(1,j), -Q1(2,j)^2-Q1(3,j)^2];
    Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
        K1(1), 0, K1(3), -K1(2);
        K1(2), -K1(3), 0, K1(1);
        K1(3), K1(2), -K1(1), 0 ] * Q1(:,j)/2;
    r_X1 = R1*Gamma1;
    Q1(:,j+1)=Q1(:,j) + Q_X1*Param.dX;
    r1(:,j+1,1)=r1(:,j,1) + r_X1*Param.dX;
end

%% Projection of the stiffness and damping in the function basis
Y=0:Param.dX:Param.L;
SK=zeros(Param.n*Param.na,Param.n*Param.na);
SD=zeros(Param.n*Param.na,Param.n*Param.na);
for k=1:Param.n_X+1
    phival=Phi(Param.na,Param.n,Y(k),Param.L); %function basis
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

%% Time integration
Param.Y=0:Param.dX:Param.L;
Kp=0.001;
Kd=0.001;
y_ref=-0.05;
yd_ref=0;
for i=1:Param.n_t %Time updating
    %Param.tau1=-0*sin(40*i*Param.dt);%-0.016256450483610;
    qa0=qa_d6(:,i);
    disp (i*Param.dt)
    init=[qa(i);q1(:,i)];
%      options=optimset('fsolve');
%      options.MaxFunEvals =1000000000000;
%      options.MaxIter = 50000000000;
%      options.TolFun=1e-5;
%      options.Jacobian='on';
% %     % options.StepTolerance=10^-4;
% %     %options.exitflag=4;
% %     %options.OptimalityTolerance=10^(-2);
%      options.Display = 'iter';
%     %options.Algorithm='levenberg-marquardt';

    options  =  optimoptions( 'fsolve', ...
                          'FiniteDifferenceStepSize', 1e-5, ...
                          'OptimalityTolerance', 1e-5, ...
                          'StepTolerance', 1e-5, ...
                          'FunctionTolerance', 1e-5, ...
                          'SpecifyObjectiveGradient', true, ...
                          'Display', 'iter',...
                          'Algorithm','levenberg-marquardt');
   global G A B C w
   C;
   G;
   w;
    tic
    %Param.Forces_Tendons=
    [x,f,a,b,J2]=fsolve(@(var)Dynamic(var,q1(:,i),qd1(:,i),qa0,eta0(:,i),i,Param),init,options);%2nd derivative of q 
    %Param.tau1=0;
    toc
    %Param.tau1=0;
    qa(:,i+1)=x(1);
    q1(:,i+1)=x(2:Param.n*Param.na+1);
    qd1(:,i+1)=(q1(:,i+1)-q1(:,i))/Param.dt;
    qa_d6(3,i+1)=qa(:,i+1);
    eta0(:,i+1)=(qa_d6(:,i+1)-qa_d6(:,i))/Param.dt;
    %init=x;
    Q1=zeros(4,Param.n_X+1);
    Q1(:,1)=[cos(qa(1,i+1)/2);0;0;sin(qa(1,i+1)/2)];
    for j=1:Param.n_X
        phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
        xia1=Param.xi_a0 + phi*q1(:,i+1);
        xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
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
        r1(:,j+1,i+1)=r1(:,j,i+1) + r_X1*Param.dX;
    end
% plot(r1(1,:,i),r1(2,:,i),'b','LineWidth', 4);axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
%     xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
%     hold on
% p1 = nsidedpoly(1000, 'Center', [r1(1,1),r1(2,1)], 'Radius', 0.01);
% plot(p1, 'FaceColor', 'g') 
% hold off
e=r1(2,end,i+1)-y_ref
r1d=(r1(2,end,i+1)-r1(2,end,i))/Param.dt;
%ed=(e-e_p)/Param.dt;
ed=r1d;
err=Kp*e+Kd*ed;
Param.tau1=(Param.tau1-err);
end
toc
%% Results view
%Robot behaviour
%for j=1:Param.n_t+1
%    j
    
%end

%End effector positions
%l_ana=zeros(Param.n_t+1,1);
m_ana=zeros(Param.n_t+1,1);
%n_ana=zeros(Param.n_t+1,1);
t=(0:Param.dt:Param.duree)';
for j=1:Param.n_t+1
    l_ana(j,1)=r1(1,end,j);
    m_ana(j,1)=r1(2,end,j);
    n_ana(j,1)=r1(3,end,j);
end
figure; hold on
grid
set(gca,'FontSize',14);
%plot(t',l_ana,'r')
plot(t',m_ana,'b')
%plot(t',n_ana,'g')
ylabel('$y$ [m] ','interpreter', 'latex');
xlabel('$t$ [s]','interpreter', 'latex');
%legend({'$x_{eff}$', '$y_{eff}$', '$z_{eff}$'}, 'interpreter', 'latex')
%title('End effector positions')

j=1;
 plot(r1(1,:,j),r1(2,:,j),'b','LineWidth', 4);axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
  p1 = nsidedpoly(1000, 'Center', [r1(1,1),r1(2,1)], 'Radius', 0.01);
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    hold on
    plot(p1, 'FaceColor', 'g') 
% gifim=[]
im={}; 
[im{j},map]=frame2im(getframe);
for j=2:Param.n_t+1
    figure(2)
    clf
    disp(j);
   plot(r1(1,:,j),r1(2,:,j),'b','LineWidth', 4);axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
   hold on
    plot(p1, 'FaceColor', 'g') 
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    [im{j},map]=frame2im(getframe);
hold off;
 pause(0.01)
end
[temp,map]=rgb2ind(im{1},Param.n_t+1);
for j=1:Param.n_t+1
    gifim(:,:,1,j)=rgb2ind(im{j},map);
end
imwrite(gifim,map,'test4_closedloop.gif');



%-------- INITIALISATION --------
figure(2); 
set(gcf,'Color','w');
set(gcf,'Position',[100 100 900 700]);

indices = [20 40 60 80 size(r1,2)];

gifname = 'Closed_loop.gif';  % <--- unified name

% ---------- PREMIERE FRAME ----------
clf; hold on;

j = 1;

plot(r1(1,:,j), r1(2,:,j),'b','LineWidth',4);

axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);
grid on; daspect([1 1 1]);
set(gca,'FontSize',14);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

% ---- Corrected indexing for poly center ----
p1 = nsidedpoly(1000, 'Center', [r1(1,1,j), r1(2,1,j)], 'Radius', 0.01);
plot(p1, 'FaceColor', 'g')

F = getframe(gcf);
[im,map] = rgb2ind(frame2im(F),256);

imwrite(im,map,gifname,'gif','LoopCount',inf,'DelayTime',0.01);



% ---------- BOUCLE PRINCIPALE ----------
for j = 2 : Param.n_t+1

    clf; hold on;

    plot(r1(1,:,j), r1(2,:,j),'b','LineWidth',4);

    axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);
    grid on; daspect([1 1 1]);
    set(gca,'FontSize',18);
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

    % ---- Corrected indexing again ----
    p1 = nsidedpoly(1000, 'Center', [r1(1,1,j), r1(2,1,j)], 'Radius', 0.01);
    plot(p1, 'FaceColor', 'g');

    F = getframe(gcf);
    im = rgb2ind(frame2im(F), map);

    imwrite(im,map,gifname,'gif','WriteMode','append','DelayTime',1);
end
