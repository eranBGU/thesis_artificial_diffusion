function [cmapmat,noObs_idx] = initialConcentrationMap(T,x,y,in_obs,c_bound,c_init)
% This function find no-obsticles indexes and the inital concentration map.

nT = size(T.loc,1);
nx = size(x,2);
ny = size(y,2);
% Calculate all the indexes for conc calculation:
noObs_idx = ones(nx,ny,nT);  % Remove Obstacles/Borders/Targets from the calculation node list
noObs_idx(1:length(x),1,:) = 0;                % borders
noObs_idx(1:length(x),end,:) = 0;              % borders
noObs_idx(1,1:length(y),:) = 0;                % borders
noObs_idx(end,1:length(y),:) = 0;              % borders
noObs_idx(repmat(in_obs,nT,1))=0;              % Obstacles

% initial Concentration Map:
cmapmat = c_bound.*ones(numel(x),numel(y),nT); % initial Concentration Map
for imap = 1:nT
    cmapmat(T.locidx(imap,1),T.locidx(imap,2),imap)=T.conc;
    noObs_idx(T.locidx(imap,1),T.locidx(imap,2),imap) = 0;  % target location
end
noObs_idx = find(noObs_idx ~= 0);         %finding all the indexs without boundary condition
cmapmat(noObs_idx) = c_init;
end

