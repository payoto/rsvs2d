

function [] = include_SnakeParam()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%% From Snakes

function [unstructReshape]=ModifUnstructured(unstructured)
    % Reshapes the unstructureddata structure to b ein line with the shape
    % of "snakes"
    unstrucFields=fieldnames(unstructured);
    nFields=length(unstrucFields);
    
    for ii=1:nFields
        field1Fields=fieldnames(unstructured.(unstrucFields{ii}));
        nFields1=length(field1Fields);
        nObjects=length(unstructured.(unstrucFields{ii}).index);
        
        for jj=1:nObjects
            for kk=1:nFields1
                if ~isstruct(unstructured.(unstrucFields{ii}).(field1Fields{kk}))
                    
                    unstructReshape.(unstrucFields{ii})(jj).(field1Fields{kk})=...
                        unstructured.(unstrucFields{ii}).(field1Fields{kk})(jj,:);
                else
                    field2Fields=fieldnames(unstructured.(unstrucFields{ii}).(field1Fields{kk}));
                    nFields2=length(field2Fields);
                    
                    for ll=1:nFields2
                        unstructReshape.(unstrucFields{ii})(jj).(...
                            [field1Fields{kk},field2Fields{ll}])=...
                            unstructured.(unstrucFields{ii}).(...
                            field1Fields{kk}).(field2Fields{ll})(jj,:);
                    end
                end
            end
        end
    end
end

function [unstructured]=ModifReshape(unstructReshape)
    % Reshapes the unstructureddata structure to b ein line with the shape
    % of "snakes"
    unstrucFields=fieldnames(unstructReshape);
    nFields=length(unstrucFields);
    
    for ii=1:nFields
        field1Fields=fieldnames(unstructReshape.(unstrucFields{ii}));
        nFields1=length(field1Fields);
        nObjects=length(unstructReshape.(unstrucFields{ii}));
        
        
        for kk=1:nFields1
            unstructured.(unstrucFields{ii}).(field1Fields{kk})=...
                zeros([nObjects,length(unstructReshape.(unstrucFields{ii})(1).(field1Fields{kk}))]);
            for jj=1:nObjects
                %                 if ~isstruct(unstructReshape.(unstrucFields{ii}).(field1Fields{kk}))
                
                unstructured.(unstrucFields{ii}).(field1Fields{kk})(jj,:)...
                    =unstructReshape.(unstrucFields{ii})(jj).(field1Fields{kk});
                %                 else
                %                     field2Fields=fieldnames(unstructReshape.(unstrucFields{ii}).(field1Fields{kk}));
                %                     nFields2=length(field2Fields);
                %
                %                     for ll=1:nFields2
                %
                %                         unstructReshape.(unstrucFields{ii}).(...
                %                             field1Fields{kk}).(field2Fields{ll})(jj,:)= ...
                %                             unstructured.(unstrucFields{ii})(jj).(...
                %                             [field1Fields{kk},field2Fields{ll}]);
                %                     end
                %                 end
            end
        end
    end
end

function [leftMost]=LeftMostCorner(coord)
    % Returns the left most coordinate in a a set
    
    [xMin]=min(coord(:,1));
    iXMin=find(coord(:,1)==xMin);
    [~,iYMin]=min(coord(iXMin,2));
    leftMost=iXMin(iYMin);
    
end

function [isCCW]=CCWLoop(coord)
    % Checks if the order of points at the left most corner to determine the
    % direction of the loop.
    [mCoord,~]=size(coord);
    %coord(end-1:end,:)=[];
    
    [leftMostCorner]=LeftMostCorner(coord);
    switch leftMostCorner
        case 1
            precVert=mCoord;
            nextVert=leftMostCorner+1;
        case mCoord
            precVert=leftMostCorner-1;
            nextVert=1;
        otherwise
            precVert=leftMostCorner-1;
            nextVert=leftMostCorner+1;
    end
    
    precVec=coord(precVert,:)-coord(leftMostCorner,:);
    nextVec=coord(nextVert,:)-coord(leftMostCorner,:);
    precAngle=ExtractAngle360([-1 -1],precVec);
    nextAngle=ExtractAngle360([-1 -1],nextVec);
    
    
    if precAngle>nextAngle
        isCCW=true;
    elseif precAngle<nextAngle
        isCCW=false;
    else
        isCCW=[];
    end
    
end

function cellSimilar=FindIdenticalVector(blockSegments)
    % this function takes in a group of Segments and returns the indices of
    % those identical grouped within a cell array. blockSegments should be a
    % vertical array of horizontal vectors to be compared.
    
    [m,n]=size(blockSegments);
    blockSegments=sort(blockSegments,2);
    preSortIndex=1:m; % save the index before shuffling
    % Shuffles edges such that similar edges are side by side
    
    for ii=1:n
        [blockSegments,sortIndex]=SortVecColumn(blockSegments,ii);
        preSortIndex=preSortIndex(sortIndex);
    end
    %compares neighbouring segments
    blockSegTrunc1=blockSegments(1:end-1,:);
    blockSegTrunc2=blockSegments(2:end,:);
    isPrecedent=[0;(sum(blockSegTrunc1==blockSegTrunc2,2)==n)];
    % creates a cell array wih as many elements as there are different edges
    cellSimilar{-sum((isPrecedent-1))}=[];
    kk=0;
    for ii=1:m
        if ~isPrecedent(ii)
            kk=kk+1;
            jj=0;
        end
        jj=jj+1;
        % assigns the presorted index to the similarity array
        cellSimilar{kk}(jj)=preSortIndex(ii);
    end
    
end

function cellSimilar=FindIdenticalVectorOrd(blockSegments)
    % this function takes in a group of Segments and returns the indices of
    % those identical grouped within a cell array. blockSegments should be a
    % vertical array of horizontal vectors to be compared.
    
    [m,n]=size(blockSegments);
    preSortIndex=1:m; % save the index before shuffling
    % Shuffles edges such that similar edges are side by side
    
    for ii=1:n
        [blockSegments,sortIndex]=SortVecColumn(blockSegments,ii);
        preSortIndex=preSortIndex(sortIndex);
    end
    %compares neighbouring segments
    blockSegTrunc1=blockSegments(1:end-1,:);
    blockSegTrunc2=blockSegments(2:end,:);
    isPrecedent=[0;(sum(blockSegTrunc1==blockSegTrunc2,2)==n)];
    % creates a cell array wih as many elements as there are different edges
    cellSimilar{-sum((isPrecedent-1))}=[];
    kk=0;
    for ii=1:m
        if ~isPrecedent(ii)
            kk=kk+1;
            jj=0;
        end
        jj=jj+1;
        % assigns the presorted index to the similarity array
        cellSimilar{kk}(jj)=preSortIndex(ii);
    end
    
end
% Order BlockEdges might cause problems as the different versions were not
% consistant.
function [cellOrderedVertex,cellOrderedEdges]=...
        OrderBlockEdges(blockEdges,blockCellTrunc)
    
    if ~exist('blockCellTrunc','var')
        blockCellTrunc=ones([length(blockEdges(:,1)),1]);
    end
    
    [mBE,~]=size(blockEdges);
    blockEdgesWorking=blockEdges;
    blockCellTruncWorking=blockCellTrunc;
    edgeList=1:mBE;
    
    % New array counters
    iCell=1;
    iEdge=0;
    % Old array locations
    ii=1;
    jj=1;
    while ~isempty(blockEdgesWorking)
        iEdge=iEdge+1;
        kk=abs(jj-3); % opposite column of jj
        % Save current Edge
        currentVertex=blockEdgesWorking(ii,jj);
        nextVertex=blockEdgesWorking(ii,kk);
        cellOrderedVertex{iCell}(iEdge,1)=currentVertex;
        cellOrderedVertex{iCell}(iEdge,2)=nextVertex;
        cellOrderedEdges{iCell}(iEdge)=edgeList(ii);
        
        % Delete current edge and edgeList entry from working set
        edgeList(ii)=[];
        blockEdgesWorking(ii,:)=[];
        currBlockCell=blockCellTruncWorking(ii);
        %Increment the counter variables
        
        blockCellTruncWorking(ii)=[];
        [ii,jj]=find(blockEdgesWorking==nextVertex);
        if length(ii)>1
            nAct=find(blockCellTruncWorking(ii)==currBlockCell);
            ii=ii(nAct);
            jj=jj(nAct);
            disp('Loops neighbouring at corner')
        end
        if isempty(ii) % reset loop if ii is not found
            % restart from the first unassigned edge
            ii=1;
            jj=1;
            % Increment the loop number
            iCell=iCell+1;
            % Restart teh edge count
            iEdge=0;
        end
    end
    
end

function [vecAngles]=ExtractAnglepm180(baseVector,testVector)
    % This function calculates the angle between vectors
    
    toComplex=[1;0+1i];
    baseAngle=angle(baseVector*toComplex);
    vecAngles=angle(testVector*toComplex)-baseAngle;
    vecAngles(vecAngles>pi)=vecAngles(vecAngles>pi)-2*pi;
    vecAngles(vecAngles<-pi)=vecAngles(vecAngles<-pi)+2*pi;
    
    
end

function [vecAngles]=ExtractAngle360(baseVector,testVector)
    % This function calculates the angle between vectors
    
    toComplex=[1;0+1i];
    baseAngle=angle(baseVector*toComplex);
    vecAngles=angle(testVector*toComplex)-baseAngle;
    vecAngles(vecAngles>(2*pi))=vecAngles(vecAngles>(2*pi))-2*pi;
    vecAngles(vecAngles<0)=vecAngles(vecAngles<0)+2*pi;
    
    
end

function [cellCentredGrid]=CellCentredGrid(refinedGrid)
    % Returns cell centred information about the grid being used
    
    cellCentredGrid=refinedGrid.cell;
    edgeCellInfo=vertcat(refinedGrid.edge(:).cellindex);
    
    origEdgeIndex=[refinedGrid.edge(:).index];
    origVertexIndex=[refinedGrid.vertex(:).index];
    origCellIndex=[refinedGrid.cell(:).index];
    
    for ii=1:length(cellCentredGrid)
        edgeCellLog=sum((edgeCellInfo==cellCentredGrid(ii).index),2)>0;
        totEdges=sum(edgeCellLog);
        cellCentredGrid(ii).edge(1:totEdges)=refinedGrid.edge(edgeCellLog);
        cellVertex=[cellCentredGrid(ii).edge(:).vertexindex];
        cellVertex=RemoveIdenticalEntries(cellVertex');
        cellVertexSub=FindObjNum(refinedGrid.vertex,cellVertex,origVertexIndex);
        totVertices=length(cellVertex);
        cellCentredGrid(ii).vertex(1:totVertices)=refinedGrid.vertex(cellVertexSub);
        
    end
end

function [vertexCentredGrid]=VertexCentredGrid(refinedGrid)
    % Returns cell centred information about the grid being used
    
    vertexCentredGrid=refinedGrid.vertex;
    vertexCellInfo=vertcat(refinedGrid.edge(:).vertexindex);
    
    origEdgeIndex=[refinedGrid.edge(:).index];
    origVertexIndex=[refinedGrid.vertex(:).index];
    origCellIndex=[refinedGrid.cell(:).index];
    
    for ii=1:length(vertexCentredGrid)
        vertCellLog=sum((vertexCellInfo==vertexCentredGrid(ii).index),2)>0;
        totEdges=sum(vertCellLog);
        vertexCentredGrid(ii).edge(1:totEdges)=refinedGrid.edge(vertCellLog);
        vertexCell=[vertexCentredGrid(ii).edge(:).cellindex];
        vertexCell=RemoveIdenticalEntries(vertexCell');
        vertexCell(vertexCell==0)=[];
        cellVertexSub=FindObjNum(refinedGrid.vertex,vertexCell,origCellIndex);
        totCells=length(vertexCell);
        vertexCentredGrid(ii).cell(1:totCells)=refinedGrid.cell(cellVertexSub);
        
    end
end

function [quotient,left]=IntegerQuotient(a,b)
    % Divides a by b and gives the integer result and the leftover
    % Works best for positive numbers
    
    quotient=floor(a/b);
    left=a-(floor(a/b)*b);
end

function surrogatePoints=PointGeneration(ranges,N_surpoints)
    % Produces an array of points containing N_surpoints in each dimension
    % combining every point with every dimesion
    %   RANGES: is a D*2 matrix containing the lower and upper bounds of each
    %           variable
    %        N: is the number of graduations in each dimension
    
    [m_ranges,~]=size(ranges);
    
    for ii=1:m_ranges
        if ranges(ii,1)~= ranges(ii,2)
            X_inter(:,ii)=linspace(ranges(ii,1),ranges(ii,2),N_surpoints);
        else
            X_inter(1:N_surpoints,ii)=ranges(ii,1);
        end
        
    end
    
    % Generation of points for RBF generation
    [Dim,~]=size(ranges);
    
    X_RBF=[];
    for ii=1:Dim
        [m_X,~]=size(X_RBF);
        m_X=max([m_X,1]);
        inter=[X_RBF,X_inter(1,ii)*ones(m_X,1)];
        for jj=2:N_surpoints
            
            % for each partial point already in corners this loop combines all
            % the values of the subsequent variable
            inter=[inter;X_RBF,X_inter(jj,ii)*ones(m_X,1)];
        end
        X_RBF=inter;
    end
    
    surrogatePoints=X_RBF;
    
end

%% Loop building

function [loopsnaxel]=ExtractSnaxelLoops(snaxel,param)
    
    [loopsnaxel]=OrderSurfaceSnaxel(snaxel);
    [loopsnaxel]=FinishLoops(loopsnaxel,param);
end

function [loopsnaxel]=FinishLoops(loopsnaxel,param)
    
    varExtract={'edgeFinish','TEShrink','LEShrink','axisRatio'};
    [edgeFinish,TEShrink,LEShrink,axisRatio]=ExtractVariables(varExtract,param);
    
    for ii=1:length(loopsnaxel)
        
        loopsnaxel(ii).snaxel.coord(:,2)=loopsnaxel(ii).snaxel.coord(:,2)*axisRatio;
    end
    
    switch edgeFinish
        case 'shrink'
            for ii=1:length(loopsnaxel)
                loopsnaxel(ii).snaxel.coord=ShrinkEdges(loopsnaxel(ii).snaxel.coord,LEShrink,TEShrink);
            end
        case 'sharpen'
            for ii=1:length(loopsnaxel)
                [loopsnaxel(ii).snaxel.coord]=SharpenEdges(loopsnaxel(ii).snaxel.coord,TEShrink,LEShrink);
            end
        case 'none'
            
    end
    
end

function [loopsnaxel]=OrderSurfaceSnaxel(snaxel,snaxPositions)
    % function extracting the snaxels into their separate loops
    if nargin==1
        global unstructglobal
        snaxPositions=PositionSnakes(snaxel,unstructglobal);
    end
    nSnax=length(snaxel);
    blockSegments=zeros(2*nSnax,2);
    for ii=1:nSnax
        for jj=0:1
            blockSegments(2*ii-jj,:)=[snaxel(ii).index,snaxel(ii).connectivity(jj+1)];
        end
    end
    
    cellSimilar=FindIdenticalVector(blockSegments);
    for ii=1:length(cellSimilar)
        blockEdgeIndex(ii)=cellSimilar{ii}(1);
    end
    blockEdges=blockSegments(blockEdgeIndex,:);
    % Order edges into closed loops
    [cellOrderedVertex]=OrderBlockEdges(blockEdges);
    snaxIndex=[snaxel(:).index];
    for ii=1:length(cellOrderedVertex)
        loopsnaxel(ii).snaxel.index=[cellOrderedVertex{ii}(:,1)];
        loopIndices=FindObjNum(snaxel,loopsnaxel(ii).snaxel.index,snaxIndex);
        loopsnaxel(ii).snaxel.coord=vertcat(snaxPositions(loopIndices).coord);
        %loopsnaxel(ii).edge.index=isEdgeIndex(cellOrderedEdges{ii});
    end
    
end

function [points]=ShrinkEdges(points,LEShrink,TEShrink)
    % Function which allows to make sharp trailing edges and leading edges
    [~,iLE]=min(points(:,1));
    [~,iTE]=max(points(:,1));
    
    [m,~]=size(points);
    
    iLEm1=iLE-1;
    [iLEm1]=IndexMod(iLEm1,m);
    iLEp1=iLE+1;
    [iLEp1]=IndexMod(iLEp1,m);
    
    iTEm1=iTE-1;
    [iTEm1]=IndexMod(iTEm1,m);
    iTEp1=iTE+1;
    [iTEp1]=IndexMod(iTEp1,m);
    
    points(iLEm1,2)=points(iLEm1,2)-LEShrink;
    points(iLEp1,2)=points(iLEp1,2)+LEShrink;
    
    points(iTEm1,2)=points(iTEm1,2)+TEShrink;
    points(iTEp1,2)=points(iTEp1,2)-TEShrink;
    
    if points(iLEm1,2)<=points(iLE+1,2)
        points([iLEm1,iLE+1],:)=[];
    end
    if points(iTEm1,2)>=points(iTEp1,2)
        points([iTEm1,iTEp1],:)=[];
    end
    
end

function [points]=SharpenEdges(points,TEShrink,LEShrink)
    % Function which allows to make sharp trailing edges and leading edges
    [~,iLE]=min(points(:,1));
    [~,iTE]=max(points(:,1));
    
    [m,~]=size(points);
    
    iLEm1=iLE-1;
    [iLEm1]=IndexMod(iLEm1,m);
    iLEp1=iLE+1;
    [iLEp1]=IndexMod(iLEp1,m);
    iLEm2=iLE-2;
    [iLEm2]=IndexMod(iLEm2,m);
    iLEp2=iLE+2;
    [iLEp2]=IndexMod(iLEp2,m);
    
    iTEm1=iTE-1;
    [iTEm1]=IndexMod(iTEm1,m);
    iTEp1=iTE+1;
    [iTEp1]=IndexMod(iTEp1,m);
    iTEm2=iTE-2;
    [iTEm2]=IndexMod(iTEm2,m);
    iTEp2=iTE+2;
    [iTEp2]=IndexMod(iTEp2,m);
    
    if LEShrink
        [points(iLEm1,:)]=Align3points(points([iLE,iLEm2],:),points(iLEm1,:));
        [points(iLEp1,:)]=Align3points(points([iLE,iLEp2],:),points(iLEp1,:));
    end
    if TEShrink
        [points(iTEm1,:)]=Align3points(points([iTE,iTEm2],:),points(iTEm1,:));
        [points(iTEp1,:)]=Align3points(points([iTE,iTEp2],:),points(iTEp1,:));
    end
    
    %     if points(iLEm1,2)<=points(iLE+1,2)
    %          points([iLEm1,iLE+1],:)=[];
    %     end
    %     if points(iTEm1,2)>=points(iTEp1,2)
    %          points([iTEm1,iTEp1],:)=[];
    %     end
    
end

function [pointAlign]=Align3points(line,point)
    
    pointAlign=point;
    
    pointAlign(2)=(line(1,2)-line(2,2))/(line(1,1)-line(2,1))...
        *(point(1)-line(2,1))+line(2,2);
    
    if isnan(pointAlign(2)) || ~isfinite(pointAlign(2))
        pointAlign=point;
    end
    
    
end

function [indMod]=IndexMod(ind,m)
    
    indMod=mod(ind-1,m)+1;
    
end

function [snakposition]=PositionSnakes(snaxel,unstructured)
    % Returns an array with Snaxel coordinates preceded by snaxel indices
    vertIndex=unstructured.vertex.index;
    vertCoord=unstructured.vertex.coord;
    fromVertex=[snaxel(:).fromvertex];
    toVertex=[snaxel(:).tovertex];
    
    nSnaxel=length(snaxel);
    
    for ii=nSnaxel:-1:1
        iToVert=vertCoord(find(vertIndex==toVertex(ii)),:); %#ok<FNDSB> % extract vertex coordinates
        iFromVert=vertCoord(find(vertIndex==fromVertex(ii)),:); %#ok<FNDSB>
        
        snakposition(ii).index=snaxel(ii).index;
        snakposition(ii).coord=iFromVert+(iToVert-iFromVert)*snaxel(ii).d;
        snakposition(ii).vectornotnorm=(iToVert-iFromVert);
        snakposition(ii).vertInit=iFromVert;
        snakposition(ii).vector=(iToVert-iFromVert)/norm(iToVert-iFromVert);
    end
    
end

