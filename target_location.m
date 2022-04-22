function [T] = target_location(T,T_loc,x,y)
nT = size(T_loc,1);
nx = size(x,2);
ny = size(y,2);
[~,T.locidx(:,1)] = min(abs(x-T_loc(:,1)),[],2); % x closest location in the mesh
[~,T.locidx(:,2)] = min(abs(y-T_loc(:,2)),[],2); % y closest location in the mesh
T.loc(:,1) = x(T.locidx(:,1));
T.loc(:,2) = y(T.locidx(:,2));
T.vinx  = sub2ind([nx ny nT],T.locidx(:,1),T.locidx(:,2))'+nx*ny*(0:1:nT-1);
end

