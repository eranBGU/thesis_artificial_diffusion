function [] = robotMovementVideo(X,Y,cmap,pathidx,c_bound)
%UNTITLED Summary of this function goes here
%   cmap - 4D vector
filestr = strcat('C:\Users\nisse\Desktop\Studies\Master\thesis\other topics\Artificial Acoustic Field for path planning\code\video\onlineObstaclesRobotNavigation',date,'.avi');
v1 = VideoWriter(filestr);
v1.FrameRate = 5;
open(v1);
f4 = figure('Visible','off');
set(f4,'PaperPositionMode', 'auto')
set(f4,'Units','Normalized','Outerposition',[0 0 1 1]);
nT = size(cmap,3);
time = size(cmap,4);
for iv = 1:20:time
    iv
    hold on
    for iplot = 1 : nT
        subplot(3,ceil(nT/3),iplot);
        hold on
        set(gca,'ColorScale','log')
        colormap('white')
        contourf(X,Y,cmap(:,:,iplot,iv),[1,0.99999999999999999,0.99999999999,0.999999999,0.99999999,0.9999999, 0.999999,  0.9999, 0.9990, 0.9900, 0.9000, 0.0100,0.0010,-0.3000, -1.0000]*c_bound)
        plot(X(pathidx(1:iv,iplot)),Y(pathidx(1:iv,iplot)),'k','LineWidth',1.5);
        xlabel('x[m]'); ylabel('y[m]')
        axis equal
    end
    frame= getframe(f4) ;
    writeVideo(v1, frame);
    %     waitbar(iv/time);
end
close(v1);
end

