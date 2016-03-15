%% Test the template file

%clear templateGrid
%include_GridCheck
%include_SnakeParam


%[templateGrid.edge,templateGrid.vertex,templateGrid.cell]=GridInit_MEX;

CheckGrid(templateGrid)
cellIndList=vertcat(templateGrid.edge(:).cellindex);
cellIndList(cellIndList<0)=0;
vertIndList=vertcat(templateGrid.edge(:).vertexindex);
vertInd=[templateGrid.vertex(:).index];
for ii=1:length(templateGrid.cell),
    [X,Y]=find(templateGrid.cell(ii).index==cellIndList);
    cellVert=vertIndList(X,:);
    cellVertSub=FindObjNum([],cellVert(:),vertInd);
    coordCell(ii,1:2)=sum(vertcat(templateGrid.vertex(cellVertSub).coord))/length(cellVertSub);
    text(coordCell(ii,1),coordCell(ii,2),int2str(templateGrid.cell(ii).index));
end

%% Memory requirement
%{
vec=[1:100:10001];
n=1;
[xGrid,yGrid]=meshgrid(vec,vec);
mem=zeros(length(vec));
for mm=1:length(vec)
    for ll=1:length(vec)
        ii=xGrid(mm,ll);
        jj=yGrid(mm,ll);
        mem(mm,ll)=2*((4+8+4+ceil((4*n)/16)*16)*(ii*jj)...
                    +(4+4*2+4*2)*((ii+1)*jj+(jj+1)*ii)...
                    +(4+8*2)*((ii+1)*(jj+1)));
    end
end

figure
c=contour(xGrid,yGrid,mem/(1e9));
clabel(c);
%}