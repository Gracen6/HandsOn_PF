
%============================================= yield the gif by the figures

stepall = 65;  % total number of the figures

for i = 1 : stepall

    % picname, should be the same with the figure names
    picname = ['GrainBoundary_' num2str(i * 40) '.png'];

    % open(picname)
    A = imread(picname);
    [I, map] = rgb2ind(A, 256);        

    if i==1
        imwrite(I,map,'GrainBoundary.gif','gif', 'Loopcount',inf,'DelayTime',0.001); % set up at the first time
    else
        imwrite(I,map,'GrainBoundary.gif','gif','WriteMode','append','DelayTime',0.001);
    end;  

    close all

end


for i = 1 : stepall

    % picname, should be the same with the figure names
    picname = ['PolyCrys_' num2str(i * 40) '.png'];

    % open(picname)
    A = imread(picname);
    [I, map] = rgb2ind(A, 256);        

    if i==1
        imwrite(I,map,'PolyCrys.gif','gif', 'Loopcount',inf,'DelayTime',0.001); % set up at the first time
    else
        imwrite(I,map,'PolyCrys.gif','gif','WriteMode','append','DelayTime',0.001);
    end;  

    close all

end
