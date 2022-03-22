% 2D diffusion equation
% !git push %to push the folder into github
clc, clear, close all;

%% Simulation Parameters
Lx = 1;                             % plate width (m)
Ly = 1;                             % plate length (m)
nx = 80;                            % number of nodes in x direction
ny = nx;                            % number of nodes in y direction                                                                             מממממממממממממממ  מנהמצת‚צמ 
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[Y,X] = meshgrid(y,x);              % generate 2D   
% D = 1;                              % diffusion constant
% tolerence = 0.5;                   % tolerence for numerical simulation (0.5 Molar*)
error_ss = 1e-4;%0.01;               % tolerence for steady state section (0.01 Molar*)

%% Initial Conditions
c_init = 10000 ;                   % Initial concentration in all nodes ( the whole plate )
c_bound = c_init ;                 % Concentraiton on the boundary conditions ("Dirichlet Conditions" )
T.conc = -c_init;

% Targets location
T_loc(1,:) = [0.611 0.611];
T_loc(2,:) = [0.611 0.111]; % Second target
T_loc(3,:) = [0.111 0.111]; % 3rd target
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
[~,R.locidx(:,1)] = min(abs(x-R_loc(:,1)),[],2); % x closest location in the mesh
[~,R.locidx(:,2)] = min(abs(y-R_loc(:,2)),[],2); % y closest location in the mesh
R.loc(:,1) = x(R.locidx(:,1));
R.loc(:,2) = y(R.locidx(:,2));


% Obstacles location
O.borders = [0 0 Lx Lx ; 0 Ly Ly 0];
O(1).obs = polyshape([0.45 0.45 0.55 0.55],[0.03 0.8 0.8 0.03]);
O(2).obs = polyshape([0.05 0.3 0.3 0.05],[0.75 0.75 0.5 0.5]);
O(3).obs = polyshape([0.1 0.6 0.6 0.1],[0.9 0.9 0.83 0.83]);

in_obs = inpolygon(X(:), Y(:), O(1).obs.Vertices(:,1),O(1).obs.Vertices(:,2));
for iobsplot =2:length(O)
    in_obs = in_obs | inpolygon(X(:), Y(:), O(iobsplot).obs.Vertices(:,1),O(iobsplot).obs.Vertices(:,2));
end

% Calculate all the indexes 
calc_idx = ones(numel(x),numel(y),nT);
calc_idx(1:length(x),1,:) = 0;                % borders
calc_idx(1:length(x),end,:) = 0;              % borders
calc_idx(1,1:length(y),:) = 0;                % borders
calc_idx(end,1:length(y),:) = 0;              % borders
calc_idx(repmat(in_obs,nT,1))=0;                         % Obstacles

% Concentration Map
cmap = c_init.*ones(numel(x),numel(y),nT);
for imap = 1:nT
R_running_index = 1:1:nT;
R_running_index = R_running_index(R_running_index~=imap); % remove the ith robot from the list
cmap(T.locidx(imap,1),T.locidx(imap,2),imap)=T.conc;
calc_idx(T.locidx(imap,1),T.locidx(imap,2),imap) = 0;  % target location
calc_idx(R.locidx(imap,1),R.locidx(imap,2),imap) = 0;
end

run_idx = find(calc_idx ~= 0);              %finding all the indexs without boundary condition

%% Online concentrations calculation and Robots Navigation:
cnorm = zeros(1e+4,1);  % for measure the online error from SS
cnorm(1) = norm(cmap(:));
cmapvec =cmap(:);
R.vinx(1,:) = sub2ind(size(cmap),R.locidx(:,1),R.locidx(:,2))'; %locidx(3) is the vectorize index of the robot
R.vinx(1,2:end) = R.vinx(1,2:end)+nx*ny.*(1:nT-1); % converting to vector shape requires the addition of values to the robots that are greater than 1     
dir.vec = [0, 1, -1, nx, nx+1, nx-1, -nx, -nx+1, -nx-1]; % This vector determines all directions of motion of the robot

it = 1; 
isoutSS = 1;
isoutTarget = 1;
norm_iteration = 10;

while isoutTarget && it <= 1e+4 %for a big nx we need more it to get the SS
    cmapvec(run_idx) = 0.25*(cmapvec(run_idx-1)+cmapvec(run_idx+1)+cmapvec(run_idx-nx)+cmapvec(run_idx+nx));
    [~,dir.idx] = min([cmapvec(R(1).vinx(it,:))'; cmapvec(R(1).vinx(it,:)+1)'; cmapvec(R(1).vinx(it,:)-1)';...
        cmapvec(R(1).vinx(it,:)+nx)'; cmapvec(R(1).vinx(it,:)+nx+1)'; cmapvec(R(1).vinx(it,:)+nx-1)';...
        cmapvec(R(1).vinx(it,:)-nx)'; cmapvec(R(1).vinx(it,:)-nx+1)'; cmapvec(R(1).vinx(it,:)-nx-1)']);
    R(1).vinx(it+1,:) = R(1).vinx(it,:)+dir.vec(dir.idx);
    
    % adding robots as obstacles 
    % !!! adding a fitcher of end calculating concentration map after the robot
    % arive the its target!!!
    if R.vinx(it+1,:) == T.vinx
        isoutTarget = 0;
    end
    it = it+1;
    cmap(:,:,:,it) = reshape(cmapvec,nx,ny,nT);

end

R.vinx(:,2:end) = R.vinx(:,2:end)-nx*ny.*(1:nT-1); % converting to vector shape requires the addition of values to the robots that are greater than 1     

%% Making video of the robot movement
robotMovementVideo(X,Y,cmap,R.vinx)