% 2D diffusion equation
% !git push %to push the folder into github
clc, clear, close all;

%% Simulation Parameters
Lx = 40;                             % plate width (m)
Ly = 4;                             % plate length (m)
nx = 60;                            % number of nodes in x direction
ny = Ly/Lx*nx;                            % number of nodes in y direction  
np = nx*ny;                         % number of nodes  
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[Y,X] = meshgrid(y,x);              % generate 2D   
numberIterations = 1e+4;

%% Initial Conditions
c_init = 10000 ;                   % Initial concentration in all nodes ( the whole plate )
c_bound = 10000 ;                 % Concentraiton on the boundary conditions ("Dirichlet Conditions" )
T.conc = -c_bound;
car_size = [2.5 1.5];

% Targets location
%----------------------%
T_loc(1,:) = [1.5 3];
T_loc(2,:) = [4.5 3]; % Second robot
T_loc(3,:) = [7.5 3]; % 3rd robot
T_loc(4,:) = [1.5 1]; % 4th robot
T_loc(5,:) = [4.5 1];
T_loc(6,:) = [7.5 1];
 
nT = size(T_loc,1);
T = target_location(T,T_loc,x,y);
%----------------------%

% Robots location
%----------------------%
R_loc(1,:) = [Lx - 4.5, 3.5];
R_loc(2,:) = [Lx - 8.5, 3.5]; % Second target
R_loc(3,:) = [Lx - 12.5, 3.5]; % 3rd target
R_loc(4,:) = [Lx - 4.5, 0]; % 
R_loc(5,:) = [Lx - 8.5, 0]; % 
R_loc(6,:) = [Lx - 12.5, 0]; % 
R = robot_location(R_loc,x,y);
%----------------------%

% Robots shape
%----------------------%
R.size(1,:) = [car_size(1),car_size(2)]/2;
R.size(2,:) = [car_size(1),car_size(2)]/2;
R.size(3,:) = [car_size(1),car_size(2)]/2;
R.size(4,:) = [car_size(1),car_size(2)]/2;
R.size(5,:) = [car_size(1),car_size(2)]/2;
R.size(6,:) = [car_size(1),car_size(2)]/2;
[R,in_robot] = isinRobot(R,X,Y,nT);
%----------------------%

% Obstacles location
%----------------------%
O = polyshape();
% O(1) = polyshape([28 30 30 28],[3 3 4 4]);
% O(2) = polyshape([24 26 26 24],[0 0 1 1]);

in_obs = isinObstacle(O,X,Y);
%----------------------%

% Calculate no-obstacles indexes and Initical Concentraition Map: 
[cmapmat,noObs_idx] = initialConcentrationMap(T,x,y,in_obs,c_bound,c_init);

%% Online concentrations calculation and Robots Navigation:
cmapvec =cmapmat(:); % matrix to vector
R.vinx(1,:) = sub2ind(size(cmapmat),R.locidx(:,1),R.locidx(:,2))'; %robots location in matrix configuration to vectorize-robot location
R.vinx(1,2:end) = R.vinx(1,2:end)+nx*ny.*(1:nT-1); % converting to vector shape requires the addition of values to the robots that are greater than 1     
dir.vec = [0, 1, -1, nx, nx+1, nx-1, -nx, -nx+1, -nx-1]; % This vector determines all directions of motion of the robot

it = 1; 
isoutTarget = 1;
calcinrobot_idx = ones(numel(x),numel(y),nT); %PreAllocation
robotInTarget_idx =[]; %contain all the nodes inside robots that are in there targets

while isoutTarget && it <= numberIterations %for a big number of nx, it sould be bigger for SS-convergence
    % Add robots as obstacles:
    calcinrobot_idx(:) = 1;
    calcinrobot_idx(in_robot)=0;
    inrobot_idx = find(calcinrobot_idx == 0);
    cmapvec(setdiff(inrobot_idx,robotInTarget_idx)) = c_bound;  % Currently the code stops updating robots that have arrived. Not sure this is the right thing to do. To change it change the line to - cmapvec(setdiff(inrobot_idx)) = c_bound;
    run_idx =setdiff(noObs_idx,[inrobot_idx;robotInTarget_idx]); %run_idx contains all the nodes at which the concentration calculation will be made
    cmapvec(run_idx) = 0.25*(cmapvec(run_idx-1)+cmapvec(run_idx+1)+cmapvec(run_idx-nx)+cmapvec(run_idx+nx));  %concentration map gauss seidel method
    [~,dir.idx] = min([cmapvec(R.vinx(it,:))'; cmapvec(R.vinx(it,:)+1)'; cmapvec(R.vinx(it,:)-1)';...
        cmapvec(R.vinx(it,:)+nx)'; cmapvec(R.vinx(it,:)+nx+1)'; cmapvec(R.vinx(it,:)+nx-1)';...
        cmapvec(R.vinx(it,:)-nx)'; cmapvec(R.vinx(it,:)-nx+1)'; cmapvec(R.vinx(it,:)-nx-1)']);  
    R.vinx(it+1,:) = R.vinx(it,:)+dir.vec(dir.idx);  % robot navigation by gradient descent
    
    R.vreal = R.vinx; %R.vreal  show the real index of the ith-robot lovation
    R.vreal(:,2:end) = R.vinx(:,2:end)-nx*ny.*(1:nT-1); 

    check_robintar = R.vinx(it+1,:) == T.vinx; %if the ith-robots in its target the ith-check_robintar = 1
    if any(check_robintar) % if one of check_robintar elements are not 0
        i_robintar = find(check_robintar);
        robotInTarget_idx = (1:np)'+np*(i_robintar-1);
        robotInTarget_idx =robotInTarget_idx(:);
        
        if any(~check_robintar)==0 % if all check_robintar elements are 1 finish the while loop
            isoutTarget = 0;
        end
    end
    it = it+1;
    cmapmat(:,:,:,it) = reshape(cmapvec,nx,ny,nT);
    [R,in_robot] = isinRobot_online(R,X,Y,nT);
end

cmap = cmapmat(:,:,:,end);
R.vinx(:,2:end) = R.vinx(:,2:end)-nx*ny.*(1:nT-1); % converting to vector shape requires the addition of values to the robots that are greater than 1     

%% Plot 2D Mesh
figure(1)
subplot(2,1,2)
plotMesh(X,Y,R,T,O,in_obs,in_robot,nx,ny,nT)

%% %% Plotting robot path on cmap
figure(2)
cmapPath(X,Y,R,cmap,nT,c_bound)

%% Making video of the robot movement
robotMovementVideo(X,Y,cmapmat,R.vinx,c_bound)

%% Functions
function [] = plotMesh(X,Y,R,T,O,in_obs,in_robot,nx,ny,nT)
hold on
plot(X(in_obs),Y(in_obs),'r+') % points inside
plot(X(~in_obs),Y(~in_obs),'b+') % points outside
plot(X(in_robot(1:nx*ny)),Y(in_robot(1:nx*ny)),'w+') % point inside the robot

if nT > 1
    plot(X(in_robot(nx*ny+1:2*nx*ny)),Y(in_robot(nx*ny+1:2*nx*ny)),'w+') % point inside the robot
end

for iobs =1:length(O)
    plot(O(iobs),'FaceColor','red');
end

for ir = 1:nT
plot(R.shape(ir),'FaceColor','green');
end

plot(R.loc(:,1),R.loc(:,2),'ko','LineWidth',1) % Robot location
plot(T.loc(:,1),T.loc(:,2),'k*','LineWidth',1) % Target location
axis equal
xlabel('x[m]'); ylabel('y[m]');
set(gcf,'Color','w')
hold off
end

function [] = cmapPath(X,Y,R,cmap,nT,c_bound)
for iplot = 1 : nT
    subplot(3,ceil(nT/3),iplot)
    hold on
    set(gca,'ColorScale','log')
    colormap('white')
    contourf(X,Y,cmap(:,:,iplot),[1,0.99999999999999999,0.99999999999999999,0.99999999999,0.999999999,0.99999999,0.9999999, 0.999999,  0.9999, 0.9990, 0.9900, 0.9000, 0.0100,0.0010,-0.3000, -1.0000]*c_bound)
    plot(X(R.vinx(:,iplot)),Y(R.vinx(:,iplot)),'k','LineWidth',1.5);
    xlabel('x[m]'); ylabel('y[m]')
    axis equal
end
end