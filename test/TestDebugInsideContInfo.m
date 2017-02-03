gridRefined=refinedGrid;
insideContourInfo=restartsnake.insideContourInfo;
vertInd=[gridRefined.vertex(:).index];
figure, hold on,
for ii=find(insideContourInfo)';
    actCoord=[gridRefined.vertex(FindObjNum([],gridRefined.edge(ii).vertexindex,vertInd)).coord];
    plot(actCoord(1:2:end),actCoord(2:2:end),'b-');
end

%[snakposition]=PositionSnakes(snaxel,refinedGriduns);
%[snakposition]=SnaxelNormal2(snaxel,snakposition);
plotPoints= @(points) plot(points([1:end],1),points([1:end],2),'r*');
plotPoints(vertcat(snakposition(:).coord))