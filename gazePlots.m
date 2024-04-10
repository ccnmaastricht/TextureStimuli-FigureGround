function [meanX, meanY] = gazePlots(xGaze, yGaze, nBlocks, nTrialsPerBlock)


for i = 1:nBlocks
    subplot(2,5,i)
    for j = 1:length(xGaze)
        xGaze(xGaze==0) = nan;
        yGaze(yGaze==0) = nan;
        
%         scatter(nanmean(xGaze(:,3)),nanmean(yGaze(:,3)))
        X = xGaze(find(xGaze(:,1)==i),3);
        Y = yGaze(find(yGaze(:,1)==i),3);
        scatter(X, Y);
        hold on;
%         xlim([-0.05 0.05]);
%         ylim([-0.05 0.05]);
        
    end
    title(['Gaze Space for block: ',num2str(i)]);
    grid on;
    
 % get mean and std dev
    meanX(i) = nanmean(nanmean(X)); stdX(i) = std(nanmean(X));
    meanY(i) = nanmean(nanmean(Y)); stdY(i) = std(nanmean(Y));
%     errorbar([meanX meanY] , [stdX stdY], 'rx')
end    
