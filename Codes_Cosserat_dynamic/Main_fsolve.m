%Main
%Parameters
clc
clear
tic
format long
%% Mechanical characteristics
Param.E = 2.563*10^5; % Young modulu
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
Param.D=10^(-4)*eye(6);%Damping matrix
Param.duree=5;%simulation time
Param.dt=0.01; %time step
Param.dX=Param.L/100; %spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L);
Param.n_t=length(Param.dt:Param.dt:Param.duree);
Param.n=3;  %strain modes number
Param.na=3; %number of actuated strains

%% Tendons parameters
Param.Rb= 0.02; % Distance between a tendon and the backbone
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)];% coordinate of the tendon in a cross-section
Param.Tendons_list=zeros(3*(Param.n_X+1),4);
for i=1:Param.n_X+1
    % Parallel case
    Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)]; 
    
    % Convergent routing
    %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-i*Param.dX))]; 
    
    %Cross routing
    %Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi/2,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(pi,Param.Rb*(Param.L-2*i*Param.dX)),Param.Tendon_coordinate(3*pi/2,Param.Rb*(Param.L-2*i*Param.dX))]; 

    %Parallel truncated
    % if (i-1)*Param.dX<Param.L/2
    %     Param.Tendons_list(3*(i-1)+1:3*i,:) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2,Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2,Param.Rb)];
    % end
end
Param.Forces_Tendons=[0;0;0;0]; % tension in each tendon

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

Param.B     = [0 0 0 ;0 1 0;0 0 1;1 0 0;0 0 0;0 0 0];
Param.B_bar = [1 0 0;0 0 0;0 0 0;0 0 0;0 1 0;0 0 1];
Param.xi_0=[0;0;0;1;0;0];
Param.xi_a0=Param.B'*Param.xi_0;
Param.xi_c=Param.B_bar'*Param.xi_0;

Param.Ha=Param.B'*Param.H*Param.B;%the matrix of the reduced Hooke coefficients
Param.Da=Param.B'*Param.D*Param.B;%reduced damping matrix

%% Initialization
%q1=zeros(Param.na*Param.n,Param.n_t+1);%Generalized coefficients
 q1(:,1)=[   0.003765254620176
  -0.001968686076627
  -0.003076876505390
   6.424753811028415
 -11.205819028548238
   6.446175216284611
                   0
                   0
                   0];

qd1=zeros(Param.na*Param.n,Param.n_t+1);%Generalized velocities
r1=zeros(3,Param.n_X+1,Param.n_t+1);%Positions of the robot sections
Q1=zeros(4,Param.n_X+1); %Quaternions
Q1(:,1)=[1;0;0;0];
Param.g0=eye(4);
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
                          'Display', 'off');
    %toc
    %options.Jacobian='off';  
    tic
    [q1(:,i+1),f,a,b,J2]=fsolve(@(var)Dynamic_V2_grad(var,q1(:,i),qd1(:,i),Param),init,options);%2nd derivative of q    
    qd1(:,i+1)=(q1(:,i+1)-q1(:,i))/Param.dt;
    toc
    rm1=zeros(3,Param.n_X+1);
    rm1(1,1)=-Param.Rb;
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
        rm1(:,j+1)=r1(:,j+1,i+1)+R1*Param.Tendon_coordinate(pi/2, Param.Rb);
    end
end
toc
%% Results view
%Robot behaviour
for j=1:Param.n_t+1
    j
    plot3(r1(1,:,j),r1(2,:,j),r1(3,:,j),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');drawnow
end

% %End effector positions
% l_ana=zeros(Param.n_t+1,1);
% m_ana=zeros(Param.n_t+1,1);
% n_ana=zeros(Param.n_t+1,1);
% t=(0:Param.dt:Param.duree)';
% for j=1:Param.n_t+1
%     l_ana(j,1)=r1(1,end,j);
%     m_ana(j,1)=r1(2,end,j);
%     n_ana(j,1)=r1(3,end,j);
% end
% figure; hold on;
% grid
% plot(t', l_ana, 'r')
% plot(t', m_ana, 'b')
% plot(t', n_ana, 'g')


figure
clf
hold on
grid on
daspect([1 1 1])
axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
view(45, 25); % Vue en perspective
camproj('perspective'); % Active la perspective

% --- Ajout d’un repère X, Y, Z ---
quiver3(0, 0, 0, Param.L/2, 0, 0, 'r', 'LineWidth', 2); % Axe X en rouge
quiver3(0, 0, 0, 0, Param.L/2, 0, 'r', 'LineWidth', 2); % Axe Y en vert
quiver3(0, 0, 0, 0, 0, Param.L/2, 'r', 'LineWidth', 2); % Axe Z en bleu
text(Param.L/2, 0, 0, 'X', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
text(0, Param.L/2, 0, 'Y', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
text(0, 0, Param.L/2, 'Z', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');

im = {};

for j = 1:Param.n_t+1
    cla % efface les tracés précédents sans effacer les axes
    
    % --- Réaffiche le repère à chaque frame ---
    quiver3(0, 0, 0, Param.L/2, 0, 0, 'r', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, Param.L/2, 0, 'r', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 0, Param.L/2, 'r', 'LineWidth', 2);
    text(Param.L/2, 0, 0, 'X', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    text(0, Param.L/2, 0, 'Y', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    text(0, 0, Param.L/2, 'Z', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    
    % --- Trace la position actuelle ---
    plot3(r1(1,:,j), r1(2,:,j), r1(3,:,j), 'b', 'LineWidth', 4);
    %plot3([-0.1;rm1(1,1)],[rm1(2,1);rm1(2,1)],'g','LineWidth', 2);
    %plot3(rm1(1,:),rm1(2,:),rm1(3,:),'g','LineWidth', 2);
    
    grid on;
    axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    view(45, 25);
    camproj('perspective');

    % --- Chrono overlay 2D (fixe à l’écran) ---
    if exist('t','var') && length(t) >= j
        current_time = t(j);
    else
        current_time = (j-1) * Param.dt;
    end

    annotation('textbox', [0.75 0.85 0.2 0.1], ...
               'String', sprintf('%.2f s', current_time), ...
               'FontSize', 14, 'FontWeight', 'bold', ...
               'Color', 'r', 'EdgeColor', 'none', ...
               'BackgroundColor', 'none');

    drawnow
    [im{j}, map] = frame2im(getframe(gcf));

    delete(findall(gcf, 'Type', 'annotation'));
end

% --- Création du GIF ---
[temp, map] = rgb2ind(im{1}, Param.n_t+1);
for j = 1:Param.n_t+1
    gifim(:,:,1,j) = rgb2ind(im{j}, map);
end
imwrite(gifim, map, 'Force_follower.gif');ylabel('$x$, $y$, $z$ [m]', 'interpreter', 'latex');
xlabel('$t$ [s]', 'interpreter', 'latex');
legend({'$x_{eff}$', '$y_{eff}$', '$z_{eff}$'}, 'interpreter', 'latex')
%title('End effector positions')





