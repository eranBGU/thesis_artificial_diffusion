function [R,in_robot] = isinRobot_online(R,X,Y,nT)
% Online finding the nodes inside the robots inside the concentration- loop

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
