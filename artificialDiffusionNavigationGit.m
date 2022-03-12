% 2D diffusion equation
clc, clear, close all;

%% Simulation Parameters
Lx = 1;                             % plate width (m)
Ly = 1;                             % plate length (m)
nx = 100;                         % number of nodes in x direction
ny = nx;                            % number of nodes in y direction 
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);

[Y,X] = meshgrid(y,x);              % generate 2D mesh

D = 1;                              % diffusion constant

tolerence = 0.5;                   % tolerence for numerical simulation (0.5 Molar*)
error_ss = 1e-4;%0.01;               % tolerence for steady state section (0.01 Molar*)

%% Initial Conditions
c_init = 10000 ;                   % Initial concentration in all nodes ( the whole plate )
c_bound = c_init ;                 % Concentraiton on the boundary conditions ("Dirichlet Conditions" )
T.conc = -10000;

% Targets location
T_loc(1,:) = [0.611 0.611];
T_loc(2,:) = [0.111 0.611]; % Second target
[~,T.locidx(1)] = min(abs(x-T_loc(1,1))); % x closest location in the mesh
[~,T.locidx(2)] = min(abs(y-T_loc(1,2))); % y closest location in the mesh
T.loc(1) = x(T.locidx(1));
T.loc(2) = y(T.locidx(2));

% Robots location
R_loc = [0.311 0.411];
[~,R.locidx(1)] = min(abs(x-R_loc(1))); % x closest location in the mesh
[~,R.locidx(2)] = min(abs(y-R_loc(2))); % y closest location in the mesh
R(1).loc(1) = x(R.locidx(1));
R(1).loc(2) = y(R.locidx(2));
number_robots = 1;

% Obstacles location
O.borders = [0 0 Lx Lx ; 0 Ly Ly 0];
O(1).obs = polyshape([0.45 0.45 0.55 0.55],[0.03 0.8 0.8 0.03]);
O(2).obs = polyshape([0.05 0.3 0.3 0.05],[0.75 0.75 0.5 0.5]);
O(3).obs = polyshape([0.1 0.6 0.6 0.1],[0.9 0.9 0.83 0.83]);

in_obs = inpolygon(X(:), Y(:), O(1).obs.Vertices(:,1),O(1).obs.Vertices(:,2));
for iobsplot =2:length(O)
    in_obs = in_obs | inpolygon(X(:), Y(:), O(iobsplot).obs.Vertices(:,1),O(iobsplot).obs.Vertices(:,2));
end
% Concentration Map
cmap = c_init.*ones(numel(x),numel(y),numel(number_robots));
cmap(T.locidx(:,1),T.locidx(:,2))=T.conc;

% Calculate all the indexes 
calc_idx = ones(size(X));
calc_idx(1:length(x),1) = 0;                % borders
calc_idx(1:length(x),end) = 0;              % borders
calc_idx(1,1:length(y)) = 0;                % borders
calc_idx(end,1:length(y)) = 0;              % borders
calc_idx(T.locidx(:,1),T.locidx(:,2)) = 0;  % target location
calc_idx(in_obs)=0;                         % Obstacles
run_idx = find(calc_idx ~= 0);              %finding all the indexs without boundary condition

%% Plot 2D Mesh
figure(1)
hold on
plot(X(in_obs),Y(in_obs),'r+') % points inside
plot(X(~in_obs),Y(~in_obs),'b+') % points outside

for iobsplot =1:length(O)
    plot(O(iobsplot).obs,'FaceColor','red');
end

plot(R.loc(1),R.loc(2),'go') % Robot location
plot(T.loc(1),T.loc(2),'g*') % Target location
axis equal
xlabel('x[m]'); ylabel('y[m]');
set(gcf,'Color','w')
hold off

%% SS concentration calculation
cnorm = zeros(1e+4,1); 
cnorm(1) = norm(cmap);
isoutSS = 1;
it = 1;
norm_iteration = 10;
cmapvec =cmap(:);

while isoutSS && it <= 1e+4 %for a big nx we need more it to get the SS
cmapvec(run_idx) = 0.25*(cmapvec(run_idx-1)+cmapvec(run_idx+1)+cmapvec(run_idx-nx)+cmapvec(run_idx+nx));
    
    if mod(it,norm_iteration)==0 && it > norm_iteration
        cnorm(it) = norm(cmapvec);
        if abs(cnorm(it)-cnorm(it-norm_iteration)) <= error_ss 
            isoutSS = 0;
        end
    end
    it = it+1;
end
cmap = reshape(cmapvec,nx,ny);

%% Plot navigation map
figure(2)
surf(X,Y,cmap);
% zlim([0.99999999e+4 1e+4])
xlabel('x[m]'); ylabel('y[m]');zlabel('Conc[M]')

%% Robot navigation
isoutTarget = 1;
ri = 1; %robot location iteration
R(1).vinx = sub2ind(size(cmap),R.locidx(1),R.locidx(2)); %locidx(3) is the vectorize index of the robot
dir.vec = [0; 1; -1; nx; nx+1; nx-1; -nx; -nx+1; -nx-1];
while isoutTarget
    [~,dir.idx] = min([cmapvec(R(1).vinx(ri)); cmapvec(R(1).vinx(ri)+1); cmapvec(R(1).vinx(ri)-1);...
        cmapvec(R(1).vinx(ri)+nx); cmapvec(R(1).vinx(ri)+nx+1); cmapvec(R(1).vinx(ri)+nx-1);...
        cmapvec(R(1).vinx(ri)-nx); cmapvec(R(1).vinx(ri)-nx+1); cmapvec(R(1).vinx(ri)-nx-1)]);
    R(1).vinx(ri+1,1) = R(1).vinx(ri,1)+dir.vec(dir.idx);
    if R(1).vinx(ri+1,1) == R(1).vinx(ri,1)
        isoutTarget = 0;
    end
    ri =ri+1;
end
%% Plotting robot path on cmap
figure(3)
set(gca,'ColorScale','log')
colormap('white')
contourf(X,Y,cmap,[10000,9999.9999999,9999.99999,9999.9999,9999.999,9999.99,  9999, 9990, 9900, 9000, 100,10,-3000, -10000])
hold on     
gr=plot(X(R.vinx),Y(R.vinx),'k','LineWidth',1.5);
xlabel('x[m]'); ylabel('y[m]')
clear all
