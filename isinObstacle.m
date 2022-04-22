function in_obs = isinObstacle(O,X,Y)
%finding all the nodes inside obstacles
in_obs = inpolygon(X(:), Y(:), O(1).Vertices(:,1),O(1).Vertices(:,2));
for iobs =2:length(O)
    in_obs = in_obs | inpolygon(X(:), Y(:), O(iobs).Vertices(:,1),O(iobs).Vertices(:,2));
end

end

