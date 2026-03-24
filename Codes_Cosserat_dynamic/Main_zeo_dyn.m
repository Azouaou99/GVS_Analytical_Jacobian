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
Param.g =0*[0 ; 0 ; 0 ; 0 ; 0; -9.81]; % gravity
Param.L = 0.2; % Rod length
Param.A = pi*Param.r^2; % Cross section area
Param.J1 = pi*Param.r^4/2; % Polar inertia moment
Param.J2 = pi*Param.r^4/4; % Inertia moment
Param.J3 = pi*Param.r^4/4; % Inertia moment
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia
Param.D=0*10^(-6)*eye(6);%Damping matrix
Param.duree=3;%simulation time
Param.dt=0.001; %time step
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);
Param.n_t=length(Param.dt:Param.dt:Param.duree);
Param.n=3;  %strain modes number
Param.na=1; %number of actuated strains

%% Tendons parameters
Param.Rb= 0.02; % Distance between a tendon and the backbone
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)];% coordinate of the tendon in a cross-section
Param.Tendons_list=zeros(3*(Param.n_X+1),1);
for i=1:Param.n_X+1
    % Parallel case
    %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)]; 
    Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(pi, Param.Rb)]; 
    
    % Convergent routing
    %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.dX))]; 
    
    %Cross routing
    %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.dX))]; 

    %Parallel truncated
    % if (i-1)*Param.dX<Param.L/2
    %     Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)];
    % end
end
Param.Forces_Tendons=0; % tension in each tendon

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
%q1(:,1)=  [-0.085909920572953
%   0.000127103268433
 % -0.000022572163879]

% q1(:,1)=[6.429620303917535
%  -11.212951359843338
%    6.448899705153095
%    0.003766202295381
%   -0.001967763183432
%   -0.003079962039788];

qd1=zeros(Param.na*Param.n,Param.n_t+1);%Generalized velocities
r1=zeros(3,Param.n_X+1,Param.n_t+1);%Positions of the robot sections
Q1=zeros(4,Param.n_X+1); %Quaternions
Q1(:,1)=[1;0;0;0];
Param.g0=eye(4);
    rm1=zeros(3,Param.n_X+1,Param.n_t+1);
    rm1(2,1,:)=-Param.Rb;
    rm1(2,:,1)=-Param.Rb;
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
    rm1(:,j+1,1)=r1(:,j+1,1)+R1*Param.Tendon_coordinate(pi, Param.Rb);
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
Kp=50;
Kd=10;
y_ref=-0.05;
yd_ref=0;
for i=1:Param.n_t %Time updating
    disp (i*Param.dt)
    init=q1(:,i);
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
                          'FiniteDifferenceStepSize', 1e-7, ...
                          'OptimalityTolerance', 1e-7, ...
                          'StepTolerance', 1e-7, ...
                          'FunctionTolerance', 1e-7, ...
                          'SpecifyObjectiveGradient', true, ...
                          'Display', 'iter');
    %toc
    %options.Jacobian='off';  
    tic
    global A B C G w
    [q1(:,i+1),f,a,b,J2]=fsolve(@(var)Dynamic_zero_dyn(var,q1(:,i),qd1(:,i),Param),init,options);%2nd derivative of q   
    qd1(:,i+1)=(q1(:,i+1)-q1(:,i))/Param.dt;
    A    ;
    C;
    B;
    G
    toc
    cco=ctrb(A,B);
    con=rank(cco);
    obo=obsv(A,C);
    obs=rank(obo);

    w;
    %Param.Forces_Tendons=0;
    % rm1=zeros(3,Param.n_X+1,Param.n_t+1);
    % rm1(2,1,:)=-Param.Rb;
        %rm1(2,:,1)=-Param.Rb;
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
        rm1(:,j+1,i+1)=r1(:,j+1,i+1)+R1*Param.Tendon_coordinate(pi, Param.Rb);
    %rm2(:,j+1)=r2(:,j+1)+R2*Param.Tendon_coordinate(0, Param.Rb);
    end
    e=r1(2,end,i+1)-y_ref
    r1d=(r1(2,end,i+1)-r1(2,end,i))/Param.dt;
   %ed=(e-e_p)/Param.dt;
   ed=r1d;
   err=Kp*e+Kd*ed;
   Param.Forces_Tendons=(Param.Forces_Tendons+err);


    % plot(r1(1,:,i+1),r1(2,:,i+1),'b','LineWidth', 4);
    % title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    % hold on
    % xlabel('x (m)'); ylabel('y (m)');drawnow
   % plot([-0.1;rm1(1,1,i+1)],[rm1(2,1,i+1);rm1(2,1,i+1)],'g','LineWidth', 2);
   % plot(rm1(1,:,i+1),rm1(2,:,i+1),'g','LineWidth', 2);
   % %plot([-0.1;rm2(1,1)],[rm2(2,1);rm2(2,1)],'g','LineWidth', 2);
   % 
   % %plot([r1(1,1,i+1);rm1(1,1)],[r1(2,1,i+1);rm1(2,1)],'b','LineWidth', 2);drawnow;
   %  % plot([r2(1,1);rm2(1,1)],[r2(2,1);rm2(2,1)],'b','LineWidth', 2);drawnow;
   % 
   %  plot([r1(1,20,i+1);rm1(1,20,i+1)],[r1(2,20,i+1);rm1(2,20,i+1)],'b','LineWidth', 2);drawnow;
   %  %plot([r2(1,20);rm2(1,20)],[r2(2,20);rm2(2,20)],'b','LineWidth', 2);drawnow;
   % 
   %  plot([r1(1,40,i+1);rm1(1,40,i+1)],[r1(2,40,i+1);rm1(2,40,i+1)],'b','LineWidth', 2);drawnow;
   %  %plot([r2(1,40);rm2(1,40)],[r2(2,40);rm2(2,40)],'b','LineWidth', 2);drawnow;
   % 
   %  plot([r1(1,60,i+1);rm1(1,60,i+1)],[r1(2,60,i+1);rm1(2,60,i+1)],'b','LineWidth', 2);drawnow;
   %  %plot([r2(1,60);rm2(1,60)],[r2(2,60);rm2(2,60)],'b','LineWidth', 2);drawnow;
   % 
   %  plot([r1(1,80,i+1);rm1(1,80,i+1)],[r1(2,80,i+1);rm1(2,80,i+1)],'b','LineWidth', 2);drawnow;
   %  %plot([r2(1,80);rm2(1,80)],[r2(2,80);rm2(2,80)],'b','LineWidth', 2);drawnow;
   %  % 
   %   plot([r1(1,end,i+1);rm1(1,end,i+1)],[r1(2,end,i+1);rm1(2,end,i+1)],'b','LineWidth', 2);drawnow;
   %   %plot([r2(1,end);rm2(1,end)],[r2(2,end);rm2(2,end)],'b','LineWidth', 2);drawnow;
     hold off
end
toc
%% Results view
%Robot behaviour
% for j=1:Param.n_t+1
%     j
% 
% end
% 
% %End effector positions
%l_ana=zeros(Param.n_t+1,1);
m_ana=zeros(Param.n_t+1,1);
%n_ana=zeros(Param.n_t+1,1);
t=(0:Param.dt:Param.duree)';
for j=1:Param.n_t+1
 %   l_ana(j,1)=r1(1,end,j);
    m_ana(j,1)=r1(2,end,j);
  %  n_ana(j,1)=r1(3,end,j);
end
figure; hold on
grid
%plot(t',l_ana,'r')
plot(t',m_ana,'b')
%plot(t',n_ana,'g')
ylabel('$y$ [m] ','interpreter', 'latex');
xlabel('$t$ [s]','interpreter', 'latex');
%legend({'$x_{eff}$', '$y_{eff}$', '$z_{eff}$'}, 'interpreter', 'latex')
%title('End effector positions')

j=1;
i=j-1
 plot(r1(1,:,j),r1(2,:,j),'b','LineWidth', 4);axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    hold on
   plot([-0.1;rm1(1,1,i+1)],[rm1(2,1,i+1);rm1(2,1,i+1)],'g','LineWidth', 2);drawnow
   plot(rm1(1,:,i+1),rm1(2,:,i+1),'g','LineWidth', 2);drawnow
   %plot([-0.1;rm2(1,1)],[rm2(2,1);rm2(2,1)],'g','LineWidth', 2);

   plot([r1(1,1,i+1);rm1(1,1,i+1)],[r1(2,1,i+1);rm1(2,1,i+1)],'b','LineWidth', 2);drawnow;
    % plot([r2(1,1);rm2(1,1)],[r2(2,1);rm2(2,1)],'b','LineWidth', 2);drawnow;

    plot([r1(1,20,i+1);rm1(1,20,i+1)],[r1(2,20,i+1);rm1(2,20,i+1)],'b','LineWidth', 2);drawnow;
    %plot([r2(1,20);rm2(1,20)],[r2(2,20);rm2(2,20)],'b','LineWidth', 2);drawnow;

    plot([r1(1,40,i+1);rm1(1,40,i+1)],[r1(2,40,i+1);rm1(2,40,i+1)],'b','LineWidth', 2);drawnow;
    %plot([r2(1,40);rm2(1,40)],[r2(2,40);rm2(2,40)],'b','LineWidth', 2);drawnow;

    plot([r1(1,60,i+1);rm1(1,60,i+1)],[r1(2,60,i+1);rm1(2,60,i+1)],'b','LineWidth', 2);drawnow;
    %plot([r2(1,60);rm2(1,60)],[r2(2,60);rm2(2,60)],'b','LineWidth', 2);drawnow;

    plot([r1(1,80,i+1);rm1(1,80,i+1)],[r1(2,80,i+1);rm1(2,80,i+1)],'b','LineWidth',2);drawnow;
im={}; 
[im{j},map]=frame2im(getframe);
for j=2:Param.n_t+1
    i=j-1
    figure(2)
    clf
    disp(j);
   plot(r1(1,:,j),r1(2,:,j),'b','LineWidth', 4);axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
       hold on
   plot([-0.1;rm1(1,1,i+1)],[rm1(2,1,i+1);rm1(2,1,i+1)],'g','LineWidth', 2);
   plot(rm1(1,:,i+1),rm1(2,:,i+1),'g','LineWidth', 2);
   %plot([-0.1;rm2(1,1)],[rm2(2,1);rm2(2,1)],'g','LineWidth', 2);

   %plot([r1(1,1,i+1);rm1(1,1)],[r1(2,1,i+1);rm1(2,1)],'b','LineWidth', 2);drawnow;
    % plot([r2(1,1);rm2(1,1)],[r2(2,1);rm2(2,1)],'b','LineWidth', 2);drawnow;

    plot([r1(1,20,i+1);rm1(1,20,i+1)],[r1(2,20,i+1);rm1(2,20,i+1)],'b','LineWidth', 2);drawnow;
    %plot([r2(1,20);rm2(1,20)],[r2(2,20);rm2(2,20)],'b','LineWidth', 2);drawnow;

    plot([r1(1,40,i+1);rm1(1,40,i+1)],[r1(2,40,i+1);rm1(2,40,i+1)],'b','LineWidth', 2);drawnow;
    %plot([r2(1,40);rm2(1,40)],[r2(2,40);rm2(2,40)],'b','LineWidth', 2);drawnow;

    plot([r1(1,60,i+1);rm1(1,60,i+1)],[r1(2,60,i+1);rm1(2,60,i+1)],'b','LineWidth', 2);drawnow;
    %plot([r2(1,60);rm2(1,60)],[r2(2,60);rm2(2,60)],'b','LineWidth', 2);drawnow;

    plot([r1(1,80,i+1);rm1(1,80,i+1)],[r1(2,80,i+1);rm1(2,80,i+1)],'b','LineWidth', 2);drawnow;
    %plot([r2(1,80);rm2(1,80)],[r2(2,80);rm2(2,80)],'b','LineWidth', 2);drawnow;
    % 
     plot([r1(1,end,i+1);rm1(1,end,i+1)],[r1(2,end,i+1);rm1(2,end,i+1)],'b','LineWidth', 2);drawnow;
     %plot([r2(1,end);rm2(1,end)],[r2(2,end);rm2(2,end)],'b','LineWidth', 2);drawnow;
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
    [im{j},map]=frame2im(getframe);
hold off;
 pause(0.01)
end
[temp,map]=rgb2ind(im{1},Param.n_t+1);
for j=1:Param.n_t+1
    gifim(:,:,1,j)=rgb2ind(im{j},map);
end
imwrite(gifim,map,'open_loop.gif');

%-------- INITIALISATION --------
figure(2); 
set(gcf,'Color','w');
set(gcf,'Position',[100 100 900 700]);

indices = [20 40 60 80 size(r1,2)];

% ---------- PREMIERE FRAME ----------
clf; hold on;

j = 1; i = j-1;

plot(r1(1,:,j), r1(2,:,j),'b','LineWidth',4);
axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);
grid on; daspect([1 1 1]);
set(gca,'FontSize',14);     % <--- increase axis numbers
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

plot([-0.1; rm1(1,1,i+1)], [rm1(2,1,i+1); rm1(2,1,i+1)], 'g','LineWidth',2);
plot(rm1(1,:,i+1), rm1(2,:,i+1), 'g','LineWidth',2);

for k = indices
    plot([r1(1,k,i+1); rm1(1,k,i+1)], [r1(2,k,i+1); rm1(2,k,i+1)], 'b','LineWidth',2);
end

F = getframe(gcf);
[im,map] = rgb2ind(frame2im(F),256);

imwrite(im,map,'test_slide.gif','gif','LoopCount',inf,'DelayTime',0.01);



% ---------- BOUCLE PRINCIPALE ----------
for j = 2:Param.n_t+1
    i = j-1;

    clf; hold on;

    plot(r1(1,:,j), r1(2,:,j),'b','LineWidth',4);

    axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);
    grid on; daspect([1 1 1]);
    set(gca,'FontSize',18);     % <--- increase axis numbers
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

    plot([-0.1; rm1(1,1,i+1)], [rm1(2,1,i+1); rm1(2,1,i+1)], 'g','LineWidth',2);
    plot(rm1(1,:,i+1), rm1(2,:,i+1), 'g','LineWidth',2);

    for k = indices
        plot([r1(1,k,i+1); rm1(1,k,i+1)], [r1(2,k,i+1); rm1(2,k,i+1)], 'b','LineWidth',2);
    end

    F = getframe(gcf);
    im = rgb2ind(frame2im(F), map);

    imwrite(im,map,'test_slide.gif','gif','WriteMode','append','DelayTime',1);
end


% 
