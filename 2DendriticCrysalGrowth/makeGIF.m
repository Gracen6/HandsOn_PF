%============creat movie by reading structures===============

stepall = 100;  % total number of the figures

for i = 1 : stepall
    
    % picname, should be the same with the figure names
    picname = [num2str(i * 25) '.jpg'];
    
    % open(picname)
    A = imread(picname);
    [I, map] = rgb2ind(A, 256);
    
    if i==1
        imwrite(I,map,'DendriticCrystalGrowth.gif','gif',...
            'Loopcount',inf,'DelayTime',0.001); % set up for the 1st frame !
    else
        imwrite(I,map,'DendriticCrystalGrowth.gif','gif',...
            'WriteMode','append','DelayTime',0.001);
    end
    
    close all
    
end
