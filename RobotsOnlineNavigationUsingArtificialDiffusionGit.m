% 2D diffusion equation
% !git push %to push the folder into github
clc, clear, close all;

%% Simulation Parameters
Lx = 1;                             % plate width (m)
Ly = 1;                             % plate length (m)
nx = 60;                            % number of nodes in x direction
ny = nx;                            % number of nodes in y direction  
np = nx*ny;                         % number of nodes  
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[Y,X] = meshgrid(y,x);              % generate 2D   
% D = 1;                              % diffusion constant
% tlerence = 0.5;                   % tolerence for numerical simulation (0.5 Molar*)
error_oss = 1e-4;%0.01;               % tolerence for steady state section (0.01 Molar*)

%% Initial Conditions
c_init = 10000 ;                   % Initial concentration in all nodes ( the whole plate )
c_bound = c_init ;                 % Concentraiton on the boundary conditions ("Dirichlet Conditions" )
T.conc = -c_init;

% Targets location
T_loc(1,:) = [0.611 0.611];
T_loc(2,:) = [0.611 0.111]; % Second target
T_loc(3,:) = [0.111 0.111]; % 3rd target
% T_loc(4,:) = [0.021 0.021]; % 3rd target
nT = size(T_loc,1);
[~,T.locidx(:,1)] = min(abs(x-T_loc(:,1)),[],2); % x closest location in the mesh
[~,T.locidx(:,2)] = min(abs(y-T_loc(:,2)),[],2); % y closest location in the mesh
T.loc(:,1) = x(T.locidx(:,1));
T.loc(:,2) = y(T.locidx(:,2));
T.vinx = sub2ind([nx ny nT],T.locidx(:,1),T.locidx(:,2))'+nx*ny*(0:1:nT-1);

% Robots location
R_loc = [0.311 0.411];
R_loc(2,:) = [0.311 0.111]; % Second robot
R_loc(3,:) = [0.911 0.001]; % 3rd robot
% R_loc(4,:) = [0.901 0.901]; % 3rd robot
[~,R.locidx(:,1)] = min(abs(x-R_loc(:,1)),[],2); % x closest location in the mesh
[~,R.locidx(:,2)] = min(abs(y-R_loc(:,2)),[],2); % y closest location in the mesh
R.loc(:,1) = x(R.locidx(:,1));
R.loc(:,2) = y(R.locidx(:,2));

% Robots shape
R.size(1,:) = [0.1 0.14]/2;
R.size(2,:) = [0.1 0.12]/2;
R.size(3,:) = [0.1 0.15]/2;
% R.size(4,:) = [0.05 0.05]/2;
[R,in_robot] = isinRobot(R,X,Y,nT);


% Obstacles location
O.borders = [0 0 Lx Lx ; 0 Ly Ly 0];
O(1).obs = polyshape([0.45 0.45 0.55 0.55],[0.03 0.8 0.8 0.03]);
O(2).obs = polyshape([0.05 0.3 0.3 0.05],[0.75 0.75 0.5 0.5]);
O(3).obs = polyshape([0.1 0.6 0.6 0.1],[0.9 0.9 0.83 0.83]);

in_obs = inpolygon(X(:), Y(:), O(1).obs.Vertices(:,1),O(1).obs.Vertices(:,2));
for iobs =2:length(O)
    in_obs = in_obs | inpolygon(X(:), Y(:), O(iobs).obs.Vertices(:,1),O(iobs).obs.Vertices(:,2));
end

% Calculate all the indexes 
noObs_idx = ones(numel(x),numel(y),nT);  % Remove Obstacles/Borders/Targets from the calculation node list
noObs_idx(1:length(x),1,:) = 0;                % borders
noObs_idx(1:length(x),end,:) = 0;              % borders
noObs_idx(1,1:length(y),:) = 0;                % borders
noObs_idx(end,1:length(y),:) = 0;              % borders
noObs_idx(repmat(in_obs,nT,1))=0;                         % Obstacles
% calc_idx(in_robot)=0;

% Concentration Map
cmapmat = c_init.*ones(numel(x),numel(y),nT);
for imap = 1:nT
R_running_index = 1:nT;
R_running_index = R_running_index(R_running_index~=imap); % remove the ith robot from the list
cmapmat(T.locidx(imap,1),T.locidx(imap,2),imap)=T.conc;
noObs_idx(T.locidx(imap,1),T.locidx(imap,2),imap) = 0;  % target location
% calc_idx(R.locidx(imap,1),R.locidx(imap,2),imap) = 0;
end
noObs_idx = find(noObs_idx ~= 0);              %finding all the indexs without boundary condition

%% Online concentrations calculation and Robots Navigation:

% cnorm = zeros(1e+4,1);  % for measure the online error from SS
% cnorm(1) = norm(cmapmat(:));
% norm_iteration = 10;
% isoutSS = 1;

cmapvec =cmapmat(:); % matrix to vector
R.vinx(1,:) = sub2ind(size(cmapmat),R.locidx(:,1),R.locidx(:,2))'; %robots location in matrix configuration to vectorize-robot location
R.vinx(1,2:end) = R.vinx(1,2:end)+nx*ny.*(1:nT-1); % converting to vector shape requires the addition of values to the robots that are greater than 1     
dir.vec = [0, 1, -1, nx, nx+1, nx-1, -nx, -nx+1, -nx-1]; % This vector determines all directions of motion of the robot

it = 1; 
isoutTarget = 1;
calcinrobot_idx = ones(numel(x),numel(y),nT); %PreAllocation
robotInTarget_idx =[]; %contain all the nodes inside robots that are in there targets

while isoutTarget && it <= 1e+4 %for a big number of nx, it sould be bigger for SS-convergence
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
plotMesh(X,Y,R,T,O,in_obs,in_robot,nx,ny,nT)

%% %% Plotting robot path on cmap
figure(2)
cmapPath(X,Y,R,cmap,nT,c_init)

%% Making video of the robot movement
robotMovementVideo(X,Y,cmapmat,R.vinx,c_init)

%% Functions
function [] = plotMesh(X,Y,R,T,O,in_obs,in_robot,nx,ny,nT)
hold on
plot(X(in_obs),Y(in_obs),'r+') % points inside
plot(X(~in_obs),Y(~in_obs),'b+') % points outside
plot(X(in_robot(1:nx*ny)),Y(in_robot(1:nx*ny)),'w+') % point inside the robot
plot(X(in_robot(nx*ny+1:2*nx*ny)),Y(in_robot(nx*ny+1:2*nx*ny)),'w+') % point inside the robot

for iobs =1:length(O)
    plot(O(iobs).obs,'FaceColor','red');
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

function [] = cmapPath(X,Y,R,cmap,nT,c_init)
for iplot = 1 : nT
    subplot(2,ceil(nT/2),iplot)
    hold on
    set(gca,'ColorScale','log')
    colormap('white')
    contourf(X,Y,cmap(:,:,iplot),[1,0.99999999999999999,0.99999999999,0.999999999,0.99999999,0.9999999, 0.999999,  0.9999, 0.9990, 0.9900, 0.9000, 0.0100,0.0010,-0.3000, -1.0000]*c_init)
    plot(X(R.vinx(:,iplot)),Y(R.vinx(:,iplot)),'k','LineWidth',1.5);
    xlabel('x[m]'); ylabel('y[m]')
    axis equal
end
end

function [R,in_robot] = isinRobot_online(R,X,Y,nT)
in_robot = zeros(numel(X),nT); % if the i-th point in the mesh is inside robot location than in_robot(i)=1
Tvec = 1:nT;
for ishape = 1:nT
    R.shape(ishape) = polyshape([X(R.vreal(end,ishape))-R.size(ishape,1), X(R.vreal(end,ishape))-R.size(ishape,1), X(R.vreal(end,ishape))+R.size(ishape,1), X(R.vreal(end,ishape))+R.size(ishape,1)],...
        [Y(R.vreal(end,ishape))-R.size(ishape,2), Y(R.vreal(end,ishape))+R.size(ishape,2),Y(R.vreal(end,ishape))+R.size(ishape,2), Y(R.vreal(end,ishape))-R.size(ishape,2)]);   
    in_robot(:,Tvec(Tvec~=ishape)) = in_robot(:,Tvec(Tvec~=ishape)) | repmat(inpolygon(X(:), Y(:), R.shape(ishape).Vertices(:,1), R.shape(ishape).Vertices(:,2)),1,max(nT-1,1)); %find points in the mesh that are inside robots locations
end
% in_robot = in_robot';
in_robot = boolean(in_robot(:));
end

function [R,in_robot] = isinRobot(R,X,Y,nT)
in_robot = zeros(numel(X),nT); % if the i-th point in the mesh is inside robot location than in_robot(i)=1
Tvec = 1:nT;
for ishape = 1:nT
    R.shape(ishape) = polyshape([R.loc(ishape,1)-R.size(ishape,1), R.loc(ishape,1)-R.size(ishape,1), R.loc(ishape,1)+R.size(ishape,1), R.loc(ishape,1)+R.size(ishape,1)],...
        [R.loc(ishape,2)-R.size(ishape,2), R.loc(ishape,2)+R.size(ishape,2), R.loc(ishape,2)+R.size(ishape,2), R.loc(ishape,2)-R.size(ishape,2)]);   
    in_robot(:,Tvec(Tvec~=ishape)) = in_robot(:,Tvec(Tvec~=ishape)) | repmat(inpolygon(X(:), Y(:), R.shape(ishape).Vertices(:,1), R.shape(ishape).Vertices(:,2)),1,max(nT-1,1)); %find points in the mesh that are inside robots locations
end
% in_robot = in_robot';
in_robot = boolean(in_robot(:));
end
