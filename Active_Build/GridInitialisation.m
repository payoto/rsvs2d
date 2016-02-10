%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%            - Grid Initialisation -
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [unstructured,loop,unstructReshape]=...
        GridInitialisation(passDomBounds,passGridSteps,passPadding,...
        typDat,loadLogical,isCheckRes,boundstr)
    % Main function for the execution of the Subdivision process
    
    
    % Defining global variables
    
    global domainBounds % domainSize is the Bounds of the design space Form [Xmin,Xmax;Ymin,Ymax;...]
    domainBounds=passDomBounds;
    
    global nGridSteps % number of steps in design domain
    nGridSteps=passGridSteps;
    
    global nDim % number of dimensions
    nDim=length(domainBounds(:,1));
    
    
    global nPadding
    nPadding=passPadding;
    
%     typDat='foillo';
%     loadLogical=false;
%     isCheckRes=true;
    
    % Execution of subroutines
    if ~loadLogical
        
        [unstructured]=InitialisEdgeGrid(typDat);
        [unstructured]=EdgeProperties(unstructured);
        isEdge=unstructured.edge.(boundstr{1});
        cond=boundstr{3};
        [loop]=OrderSurfaceVertex(unstructured,isEdge,cond);
    else
        datPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.mat'];
        load(datPath);
    end

    if isCheckRes
        CheckResults(unstructured,loop)
    end
    disp('Reshaping')
    [unstructReshape]=ModifUnstructured(unstructured);
end


%% Surface Identification functions

function [loop]=OrderSurfaceVertex(unstructured,isEdge,cond)
    % function ordering the surface vertices in counter-clockwise order
    % isEdge is a logical indexing arraying informing which edges are on the
    % edge of a surface
    % cond is the optional string argument describing the condition fulfilled
    % by isEdge possible arg: '0bound', '1bound', 'intermBound'
    
    isEdgeIndex=find(isEdge);
    
    blockEdges=unstructured.edge.vertexindex(isEdgeIndex,:);
    blockCell=unstructured.edge.cellindex(isEdgeIndex,:);
    fillCell=unstructured.edge.fill(isEdgeIndex,:);
    if ~exist('cond','var'), cond='1bound';end
    switch cond
        case '0bound'
            fillCell=fillCell>0;
        case '1bound'
            fillCell=fillCell==1;
        case 'intermBound'
            fillCell(:,1)=fillCell(:,1)>fillCell(:,2);
            fillCell(:,2)=fillCell(:,2)>fillCell(:,1); 
    end
    parfor ii=1:length(fillCell(:,1))
        %colNum=find(fillCell(ii,:));
        blockCellTrunc(ii)=blockCell(ii,find(fillCell(ii,:)));
        
    end
    
    coordVertex=unstructured.vertex.coord;
    
    % Order edges into closed loops
    [cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(blockEdges,blockCellTrunc);
    
    for ii=1:length(cellOrderedVertex)
        loop(ii).vertex.index=[cellOrderedVertex{ii}(:,1);cellOrderedVertex{ii}(1:2,1)];
        loop(ii).vertex.coord=coordVertex(loop(ii).vertex.index,:);
        loop(ii).edge.index=isEdgeIndex(cellOrderedEdges{ii});
    end
    
end

function [cellOrderedVertex,cellOrderedEdges]=...
        OrderBlockEdges(blockEdges,blockCell)
    
    
    [mBE,~]=size(blockEdges);
    blockEdgesWorking=blockEdges;
    blockCellWorking=blockCell;
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
        cellOrderedCells{iCell}(iEdge)=blockCell(edgeList(ii));
        % Delete current edge and edgeList entry from working set
        edgeList(ii)=[];
        blockEdgesWorking(ii,:)=[];
        blockCellWorking(ii)=[];
        %Increment the counter variables
        
        [ii,jj]=find(blockEdgesWorking==nextVertex);
        if length(ii)>1
            oldII=cellOrderedEdges{iCell}(iEdge);
            iiIndex=cellOrderedCells{iCell}(iEdge)==blockCellWorking(ii);
            ii=ii(iiIndex);
            
            
            jj=jj(iiIndex);
            if isempty(ii)
                warning('ii is empty after cell identification this is an unlikely event in normal operations')
            end
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

function [unstructured]=EdgeProperties(unstructured)
    % Extracts fill data for each edge and classifies edges depending on
    % neighbouring cells
    
    % Extract fill information
    edgeFill=EdgeFillInformation(unstructured);
    
    % Extract indices of edges matching boundary criteria
    edgeFillSort=sort(edgeFill,2);
    
    is0=sum((edgeFillSort==0),2);
    is1=sum((edgeFillSort==1),2);
    isElse=sum((edgeFillSort~=0)&(edgeFillSort~=1),2);
    
    unstructured.edge.boundaryisHalf=(isElse==2);
    unstructured.edge.boundaryis0=(is0==1);
    unstructured.edge.boundaryis1=(is1==1);
    
    unstructured.edge.solidisIn0=(is0==2);
    unstructured.edge.solidnotIn0=(is0~=2);
    unstructured.edge.solidisIn1=(is1==2);
    
    unstructured.edge.fill=edgeFill;
%     unstructured.edge.boundary=boundary;
%     unstructured.edge.solid=solid;
    
end

function [edgeFill]=EdgeFillInformation(unstructured)
    % returns the fill information of neighbouring cells
    
    
    % decomposing structure
    fillCellDat=[0;unstructured.cell.fill];
    edgeCellIndex=unstructured.edge.cellindex+1;
    
    % detecting array size
    [mCI,nCI]=size(edgeCellIndex);
    
    % Preallocating
    edgeFill=zeros(mCI,nCI);
    
    for ii=1:nCI
        edgeFill(:,ii)=fillCellDat(edgeCellIndex(:,ii));
    end
    
end

%% Plot Functions
function []=CheckResults(unstructured,loop)
    global nDim domainBounds
    
    if nDim==2
        figh=figure;
        axh=axes;
        hold on
        
        colString='bgcmyk';
        
        isEdgeIndex=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeIndex)
            PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'bo')
        end
        
        isEdgeIndex=find(unstructured.edge.boundaryis0);
        for ii=1:length(isEdgeIndex)
            PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'b-')
        end
        
        
        isCellFull=find(unstructured.cell.fill);
        for ii=1:length( isCellFull)
            %PlotCell(figh,axh,unstructured, isCellFull(ii),'bs')
        end
        
        for ii=1:length(loop)
            [~,colIndex]=IntegerQuotient(ii,length(colString));
            colIndex=colIndex+1;
            PlotLoop(figh,axh,loop,ii,[colString(colIndex),'--'])
%            PlotSubDiv(figh,axh,loop,ii,[colString(colIndex),'-'])
        end
        
        axis equal
        axis([domainBounds(1,1:2) domainBounds(2,1:2)])
    end
    
end

function []=PlotEdge(figh,axh,unstructured,indexEdge,format)
    figure(figh)
    %axes(axh)
    
    vertices=unstructured.edge.vertexindex(indexEdge,:);
    coord=unstructured.vertex.coord(vertices,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotLoop(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).vertex.coord;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotCell(figh,axh,unstructured,indexCell,format)
    figure(figh)
    axes(axh)
    
    
    coord=unstructured.cell.coord(indexCell,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end


%% Initialisation Functions

function [unstructured]=InitialisEdgeGrid(typDat)
    % Main function for the execution of the Subdivision process
    
    [unstructured]=Initialisation_Square(typDat);
    edgeTemplate=EdgeBuildTemplate(unstructured);
    unstructured.edge=CreateEdges(unstructured,edgeTemplate);
    
end

function [fill]=InputData(typDat)
    global nDim nGridSteps nPadding
    % Creates geometry information data in a cell centred array
    arraySize=ones(1,nDim)*nGridSteps;
    
    switch typDat
        case 'rand'
            fill=randi(2,arraySize)-1;
            [fill]=AddZeroLayer(fill,nPadding);
            sizeIm=size(fill);
            nGridSteps=sizeIm(1);
        case 'cyli'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'wedg'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'foil'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.jpg'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'foillo'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'foillosubdiv'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uni'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uniLo1'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uniLo2'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uniLo3'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'square'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'low5shapek'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'k',nPadding);
        case 'low5shape'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'vlofoil'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'doubleBody'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'k',nPadding);
        otherwise
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
    end
    
end

function [unstructured]=Initialisation_Square(typDat)
    % Initialise a nGridSteps*nGridSteps dataset
    
    global nDim nGridSteps domainBounds
    parametrisation.fill=InputData(typDat);
    arraySize=ones(1,nDim)*nGridSteps;
    
    if arraySize~=size(parametrisation.fill)
        error('Invalid fill array has been supplied, Review array dimensions')
    end
    
    % Finds Cell centres position
    vecNum=((0:(nGridSteps-1))+0.5)/(nGridSteps);
    vecCoord{nDim}=0;
    vertexNum=linspace(0,1,nGridSteps+1);
    vertexCoord{nDim}=0;
    for ii=1:nDim
        vecCoord{ii}=vecNum*(domainBounds(ii,2)-domainBounds(ii,1))...
            +domainBounds(ii,1);
        vertexCoord{ii}=vertexNum*(domainBounds(ii,2)-domainBounds(ii,1))...
            +domainBounds(ii,1);
    end
    [parametrisation.positionCell{1:nDim}]=ndgrid(vecCoord{1:nDim});
    [parametrisation.positionVertex{1:nDim}]=ndgrid(vertexCoord{1:nDim});
    
    fullCellsIndex=find(parametrisation.fill);
    parametrisation.fullCells=zeros(length(fullCellsIndex),nDim);
    for ii=1:nDim
        parametrisation.fullCells(:,ii)=parametrisation.positionCell{ii}(fullCellsIndex);
    end
    
    % Proper edge base parametrisation
    unstructured.cell.index=(1:(nGridSteps^nDim))';
    unstructured.vertex.index=(1:((nGridSteps+1)^nDim))';
    unstructured.cell.fill=parametrisation.fill(:);
    for ii=1:nDim
        unstructured.cell.coord(:,ii)=parametrisation.positionCell{ii}(:);
        unstructured.vertex.coord(:,ii)=parametrisation.positionVertex{ii}(:);
    end
    
    [unstructured]=CellVertexIndex(unstructured);
end

function [unstructured]=CellVertexIndex(unstructured)
    % connect information from cell to vertices
    
    global nDim
    
    ranges=[zeros(nDim,1),ones(nDim,1)];
    movVecs=PointGeneration(ranges,2);
    
    cellIndex=unstructured.cell.index;
    sizCellDat=ones(1,nDim)*length(cellIndex)^(1/nDim);
    arraySub=ind2subWrap(sizCellDat,cellIndex);
    indexCellVertex=zeros(length(cellIndex),2^nDim);
    parfor ii=1:length(cellIndex)
        
        arrayCellVertSub=ones(2^nDim,1)*arraySub(ii,:)+movVecs;
        cellCellVertSub=mat2cell(arrayCellVertSub,2^nDim,ones(1,nDim));
        indexCellVertexInterm(ii,:)=sub2ind(sizCellDat+1,cellCellVertSub{:})';
    end
    
    indexCellVertex=indexCellVertexInterm(cellIndex,:);
    unstructured.cell.vertexindex=indexCellVertex;
end

function [arraySub]=ind2subWrap(arraySize,workInd)
    % wrapper function such that ind2sub does what I need it for
    
    global nDim
    
    subCell{nDim}=0;
    [subCell{:}]=ind2sub(arraySize,workInd);
    arraySub=cell2mat(subCell);
    
end

function [segmentColumns]=EdgeBuildTemplate(unstructured)
    % builds an edge building template for each cell depending on the number of
    % dimensions
    
    global nDim
    
    % Define edge building template
    possibleSegments=zeros(sum(1:((2^nDim)-1)),2);
    jj=1;
    
    for ii=1:2^nDim-1
        possibleSegments(jj:(jj+(2^nDim-(ii))-1),:)=...
            [ones(2^nDim-(ii),1)*ii,(ii+1:2^nDim)'];
        jj=jj+2^nDim-(ii);
    end
    
    % Extract vertex data for 1 cell for analysis
    templateCellVertex=unstructured.cell.vertexindex(1,:);
    vertCoord=unstructured.vertex.coord(templateCellVertex,:);
    
    
    % Calculate middle of segment (average of two coordinates)
    centroid=mean(vertCoord);
    centroidVector=zeros(length(possibleSegments(:,1)),nDim);
    normalVector=zeros(length(possibleSegments(:,1)),nDim);
    tanVector=zeros(length(possibleSegments(:,1)),nDim);
    midPoints=zeros(length(possibleSegments(:,1)),nDim);
    isSegment=zeros(1,length(possibleSegments(:,1)));
    
    for ii=1:length(possibleSegments(:,1))
        midPoints(ii,:)=(vertCoord(possibleSegments(ii,1),:)+vertCoord(possibleSegments(ii,2),:))/2;
        centroidVector(ii,:)=midPoints(ii,:)-centroid;
        % Calculate normal vectors (only works in 2D)
        tanVector(ii,:)=(vertCoord(possibleSegments(ii,1),:)-vertCoord(possibleSegments(ii,2),:));
        normalVector(ii,:)=CalcNormVec2D(tanVector(ii,:));
        
        % check wether all points are in same half plane
        
        checkSide=(vertCoord*normalVector(ii,:)')<=(midPoints(ii,:)*normalVector(ii,:)');
        isSegment(ii)=(sum(checkSide)==(length(vertCoord(:,1))))|(sum(checkSide)==2);
        
    end
    
    % Extract valid column combinations
    segIndex=find(isSegment);
    segmentColumns=possibleSegments(segIndex,:);
    
end

function normalVector=CalcNormVec2D(tanVector)
    % Calculates a vector normal to another in 2D
    
    if tanVector(2)~=0
        normalVector(1)=1;
        normalVector(2)=-normalVector(1)*tanVector(1)/tanVector(2);
    elseif tanVector(1)~=0
        normalVector(2)=1;
        normalVector(1)=-normalVector(2)*tanVector(2)/tanVector(1);
    end
    normalVector=normalVector/norm(normalVector);
    
end

function [edge]=CreateEdges(unstructured,edgeTemplate)
    % Constructs the edges from the vertex and cell data combined with the
    % edge template.
    
    % Breaking down structure
    indexCellVertex=unstructured.cell.vertexindex;
    indexCell=unstructured.cell.index;
    % Saving Lengths of arrays
    
    [mIC,nIC]=size(indexCellVertex);
    [mET,nET]=size(edgeTemplate);
    
    segmentInfo=zeros(mET*mIC,nET+1);
    for ii=1:mET
        iStart=(ii-1)*mIC+1;
        iEnd=(ii)*mIC;
        segmentInfo(iStart:iEnd,:)=[indexCellVertex(:,edgeTemplate(ii,:)),indexCell];
    end
    
    cellSimilar=FindIdenticalVector(segmentInfo(:,1:end-1));
    indexEdge=zeros(length(cellSimilar),1);
    cellIndexEdge=zeros(length(cellSimilar),2);
    parfor ii=1:length(cellSimilar)
        %nCell(ii)=length(cellSimilar{ii});
        indexEdge(ii,1)=ii;
        %cellInfo=cellSimilar{ii};
        
        cellIndexEdge(ii,:)=...
            [segmentInfo(cellSimilar{ii},end),zeros(2-length(cellSimilar{ii}))];
        vertexindexEdge(ii,:)=segmentInfo(cellSimilar{ii}(1),1:end-1);
    end
    edge.index=indexEdge;
    edge.cellindex=cellIndexEdge;
    edge.vertexindex=vertexindexEdge;
    
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
    isPrecedent=[0;(sum(blockSegTrunc1==blockSegTrunc2,2)>1)];
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

%% Treatment of Input Images
function finishedImage=ImageProcess(imPath,imType,nPad,n)
    % Processes images into a valid input to the program
    % impath indicates the background colour: 'k' is black and 'w' is white
    
    global nGridSteps % number of steps in design domain
    if ~exist('n','var'); n=0; end
    
    preProcImage=PreProcImage(imPath);
    preProcImage=ProcessType(imType,preProcImage);
    procImage=GradProcessor(preProcImage);
    finishedImage=ResizeImage(procImage);
    [finishedImage]=AddZeroLayer(finishedImage,nPad,n);
    sizeIm=size(finishedImage);
    
    nGridSteps=sizeIm(1);
    
end

function [preProcImage]=PreProcImage(imPath)
    % Load Image and reduce it to an averaged double array from 0 to 1
    
    preProcImage=imread(imPath);
    imClass=class(preProcImage);
    numBit=str2num(regexprep(imClass,'uint',''));
    preProcImage=mean(preProcImage,3);
    
%     lvlMax=max(preProcImage(:));
%     lvlMin=min(preProcImage(:));

    lvlMax=2^numBit-1;
    lvlMin=0;
    
    if lvlMin~=lvlMax
        preProcImage=(preProcImage-lvlMin)/(lvlMax-lvlMin);
    else
        preProcImage=zeros(size(preProcImage));
    end
    
end

function [imageDat]=ProcessType(imType,imageDat)
    % Processes the background of the image (either white or black)
    
    switch imType
        case 'k'
            imBackground=0;
        case 'w'
            imBackground=1;
    end
    
    imageDat=abs(imageDat-imBackground);
    
end

function [array]=SquareArray(array)
    % transform a rectangular array into a square by padding with zeros
    
    [m,n]=size(array);
    if n>m
        array=[array;zeros(n-m,n)];
    elseif n<m
        array=[array,zeros(m,m-n)];
    end
end

function [finishedImage]=ResizeImage(procImage)
    % Resizes the image to the right size for the execution of the program
    
    procImage=TrimZeros(procImage);
    finishedImage=SquareArray(procImage);
    
    finishedImage=finishedImage';
    finishedImage=finishedImage(:,end:-1:1);
end

function [procImage]=GradProcessor(preProcImage)
    % Process gradients to ones or zeros
    
    procImage=(preProcImage);
end

function [array]=TrimZeros(array)
    % this function trims complete lines and columns of zeros from the borders
    % of an array
    
    sizArray=size(array);
    nDim=length(sizArray);
    infoArray{nDim}=[];
    
    for iDim=1:nDim
        dims=1:nDim;
        dims(dims==iDim)=[];
        workingArray=array;
        for ii=dims
            workingArray=sum(workingArray,ii);
        end
        infoArray{iDim}=squeeze(workingArray);
    end
    
    for iDim=1:nDim
        
        indexCell=[];
        indexCell{nDim}=[];
        for ii=1:nDim
            indexCell{ii}=':';
        end
        
        candidates=find(infoArray{iDim}==0);
        if ~isempty(candidates)
            dimRmv=[];
            iRmv=1;
            kkLow=0;
            condition=candidates(kkLow+1)==(kkLow+1);
            while condition
                kkLow=kkLow+1;
                dimRmv(iRmv)=kkLow;
                iRmv=iRmv+1;
                condition=(kkLow+1)<=numel(candidates);
                if condition
                    condition=candidates(kkLow+1)==(kkLow+1);
                end
                
            end
            
            kkHigh=sizArray(iDim)+1;
            kk=length(candidates)+1;
            while (candidates(kk-1)==(kkHigh-1))
                kkHigh=kkHigh-1;
                kk=kk-1;
                dimRmv(iRmv)=kkHigh;
                iRmv=iRmv+1;
                if kk==1
                    break
                end
            end
            indexCell{iDim}=dimRmv;
            array(indexCell{:})=[];
        end
    end
    
end

function [array]=AddZeroLayer(array,nPad,n)
    % Adds n layer of zeros to every side of the array
    
    sizArray=size(array);
    sizNew=sizArray+2*nPad;
    
    nDim=length(sizArray);
    for ii=1:nDim
        sizArray=size(array);
        sizePads=sizArray;
        sizePads(ii)=2*nPad+sizePads(ii);
        newArray=zeros(sizePads)+n;
        for jj=1:nDim
            indexCell{jj}=':';
        end
        indexCell{ii}=nPad+1:sizePads(ii)-nPad;
        newArray(indexCell{:})=array;
        array=newArray;
    end
end

%% General Utility Functions
function [vec,iRows]=SortVecColumn(vec,iCol)
    % Sorts according to a columns
    [~,iRows]=sort(vec(:,iCol));
    vec=vec(iRows,:);
    
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

function [quotient,left]=IntegerQuotient(a,b)
    % Divides a by b and gives the integer result and the leftover
    % Works best for positive numbers
    
    quotient=floor(a/b);
    left=a-(floor(a/b)*b);
end

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

function []=template()
    
end


%% Do not Remove function development in progress
function normVector=CalcNormalNDOptim()
    % Doesn't work In Progress
    centroid=mean(vertCoord);
    centroidVector=zeros(length(possibleSegments(:,1)),nDim);
    normalVector=zeros(length(possibleSegments(:,1)),nDim);
    tanVector=zeros(length(possibleSegments(:,1)),nDim);
    midPoints=zeros(length(possibleSegments(:,1)),nDim);
    
    for ii=1:length(possibleSegments(:,1))
        midPoints(ii,:)=(vertCoord(possibleSegments(ii,1),:)+vertCoord(possibleSegments(ii,2),:))/2;
        centroidVector(ii,:)=midPoints(ii,:)-centroid;
        
        % Optimisation parameters
        tanVector(ii,:)=(vertCoord(possibleSegments(ii,1),:)-vertCoord(possibleSegments(ii,2),:));
        % tanVector(ii,:)=tanVector(ii,:)/tanVector(ii,1);
        if norm(tanVector(ii,:))~=0
            tanVector(ii,:)=tanVector(ii,:)/norm(tanVector(ii,:));
        end
        [cij,cji]=meshgrid(centroidVector(ii,:),centroidVector(ii,:));
        
        H=(cij.*cji).*(eye(nDim)+1);
        A=-ones(1,nDim);
        
        b=-1e-8;
        intermVec=quadprog(-H,[],A,b,tanVector(ii,:),0,zeros(1,nDim),ones(1,nDim));
        if isempty(intermVec)
            disp('in if loop')
            A=ones(1,nDim);
            A(1)=-A(1);
            b=1e-8;
            intermVec=quadprog(-H,zeros(1,nDim),A,b,tanVector(ii,:),0,zeros(1,nDim),ones(1,nDim));
        end
        intermVec
        normalVector(ii,:)=intermVec;
    end
    
end
% Do not Remove function development in progress
