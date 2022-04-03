function [] = robotMovementVideo(X,Y,cmap,pathidx)
%UNTITLED Summary of this function goes here
%   cmap - 4D vector 
filestr = strcat('C:\Users\nisse\Desktop\Studies\Master\thesis\other topics\Artificial Acoustic Field for path planning\code\video\onlineObstaclesRobotNavigation',date,'.avi');
v1 = VideoWriter(filestr);
v1.FrameRate = 5;
open(v1);
f4 = figure('Visible','off');
nT = size(cmap,3);
time = size(cmap,4);
for iv = 1:4:time
    iv
    hold on
    for iplot = 1 : nT
    subplot(1,nT,iplot)
    hold on
    set(gca,'ColorScale','log')
    colormap('white')
    contourf(X,Y,cmap(:,:,iplot,iv),[10000,9999.9999999,9999.99999,9999.9999,9999.999,9999.99,  9999, 9990, 9900, 9000, 100,10,-3000, -10000]*100)
    plot(X(pathidx(1:iv,iplot)),Y(pathidx(1:iv,iplot)),'k','LineWidth',1.5);
    xlabel('x[m]'); ylabel('y[m]')
    axis equal
    end
    frame= getframe(f4) ;
    writeVideo(v1, frame);
    waitbar(iv/time);
end
close(v1);
end

