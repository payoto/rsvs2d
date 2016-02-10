

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

% Order BlockEdges might cause problems as the different versions were not
% consistant.
function [cellOrderedVertex,cellOrderedEdges]=...
        OrderBlockEdges(blockEdges,discardInput)
    
    
    [mBE,~]=size(blockEdges);
    blockEdgesWorking=blockEdges;
    
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
        
        %Increment the counter variables
        
        [ii,jj]=find(blockEdgesWorking==nextVertex);
        if length(ii)>1
            warning('ii is empty after cell identification this is an unlikely event in normal operations')
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

function [vec,iRows]=SortVecColumn(vec,iCol)
    % Sorts according to a columns
    [~,iRows]=sort(vec(:,iCol));
    vec=vec(iRows,:);
    
end

function [vectorEntries]=RemoveIdenticalEntries(vectorEntries)
    % Function which removes identical entries in a column vector
    vectorEntriesUnsort=vectorEntries;
    [vectorEntries,vectorIndex]=sort(vectorEntries);
    kk=1;
    rmvDI=[];
    for ii=2:length(vectorEntries)
        if vectorEntries(ii)==vectorEntries(ii-1)
            rmvDI(kk)=ii;
            kk=kk+1;
        end
    end
    %vectorEntries(rmvDI)=[];
    vectorIndex(rmvDI)=[];
    vectorEntries=vectorEntriesUnsort(vectorIndex);
end

function sub=FindObjNum(object,index,objInd)
    % finds the array index from a snaxel number
    if ~exist('objInd','var')
        objInd=[object(:).index];
    end
    sub=zeros(length(index),1);
    additionalSlots=0;
    for ii=1:length(index)
        
        snaxLog=objInd==index(ii);
        jj=ii+additionalSlots;
        subInter=find(snaxLog);
        if isempty(subInter)
            sub(jj)=0;
        elseif numel(subInter)>1
            sub(jj:jj+length(subInter)-1)=subInter;
            additionalSlots=additionalSlots+numel(subInter)-1;
        else
            sub(jj)=subInter;
            
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

%% From Main
%{
function [isCCW]=CCWLoop(coord)
    % Checks if the order of points at the left most corner to determine the
    % direction of the loop.
    coord(end-1:end,:)=[];
    [mCoord,~]=size(coord);
    
    [xMin]=min(coord(:,1));
    iXMin=find(coord(:,1)==xMin);
    [~,iYMin]=min(coord(iXMin,2));
    leftMostCorner=iXMin(iYMin);
    
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
    
    if coord(precVert,2)>coord(nextVert,2)
        isCCW=true;
    elseif coord(precVert,2)<coord(nextVert,2)
        isCCW=false;
    else
        isCCW=[];
    end
    
end
%}
function [quotient,left]=IntegerQuotient(a,b)
    % Divides a by b and gives the integer result and the leftover
    % Works best for positive numbers
    
    quotient=floor(a/b);
    left=a-(floor(a/b)*b);
end

%% From Grid Initialisation
%{
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
function [cellOrderedVertex,cellOrderedEdges]=...
        OrderBlockEdges2(blockEdges,blockCell)
    
    
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


%}

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
%{
function [vec,iRows]=SortVecColumn(vec,iCol)
    % Sorts according to a columns
    [~,iRows]=sort(vec(:,iCol));
    vec=vec(iRows,:);
    
end
%}

%% From GridRefinement
%{
function sub=FindObjNum(object,index,objInd)
    % finds the array index from a snaxel number
    if ~exist('objInd','var')
        objInd=[object(:).index];
    end
    sub=zeros(length(index),1);
    for ii=1:length(index)
        
        snaxLog=objInd==index(ii);
        subInter=find(snaxLog);
        if isempty(subInter)
            sub(ii)=0;
        else
            sub(ii)=subInter;
        end
    end
end

function [vectorEntries]=RemoveIdenticalEntries(vectorEntries)
    % Function which removes identical entries in a column vector
    
    [vectorEntries,vectorIndex]=sort(vectorEntries);
    kk=1;
    rmvDI=[];
    for ii=2:length(vectorEntries)
        if vectorEntries(ii)==vectorEntries(ii-1)
            rmvDI(kk)=ii;
            kk=kk+1;
        end
    end
    vectorEntries(rmvDI)=[];
    vectorIndex(rmvDI)=[];
    vectorEntries=vectorEntries(vectorIndex);
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
    
    if coord(precVert,2)>coord(nextVert,2)
        isCCW=true;
    elseif coord(precVert,2)<coord(nextVert,2)
        isCCW=false;
    else
        isCCW=[];
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
%}
%% Tecplot out 
%{
function sub=FindObjNum(object,index,objInd)
    % finds the array index from a snaxel number
    if ~exist('objInd','var')
        objInd=[object(:).index];
    end
    sub=zeros(length(index),1);
    additionalSlots=0;
    for ii=1:length(index)
        
        snaxLog=objInd==index(ii);
        jj=ii+additionalSlots;
        subInter=find(snaxLog);
        if isempty(subInter)
            sub(jj)=0;
        elseif numel(subInter)>1
            sub(jj:jj+length(subInter)-1)=subInter;
            additionalSlots=additionalSlots+numel(subInter)-1;
        else
            sub(jj)=subInter;
            
        end
    end
end

function [vectorEntries]=RemoveIdenticalEntries(vectorEntries)
    % Function which removes identical entries in a column vector
    vectorEntriesUnsort=vectorEntries;
    [vectorEntries,vectorIndex]=sort(vectorEntries);
    kk=1;
    rmvDI=[];
    for ii=2:length(vectorEntries)
        if vectorEntries(ii)==vectorEntries(ii-1)
            rmvDI(kk)=ii;
            kk=kk+1;
        end
    end
    %vectorEntries(rmvDI)=[];
    vectorIndex(rmvDI)=[];
    vectorEntries=vectorEntriesUnsort(vectorIndex);
end
%}
%% Velocities
%{
function sub=FindObjNum(object,index,objInd)
    % finds the array index from a snaxel number
    if ~exist('objInd','var')
        objInd=[object(:).index];
    end
    sub=zeros(length(index),1);
    additionalSlots=0;
    for ii=1:length(index)
        
        snaxLog=objInd==index(ii);
        jj=ii+additionalSlots;
        subInter=find(snaxLog);
        if isempty(subInter)
            sub(jj)=0;
        elseif numel(subInter)>1
            sub(jj:jj+length(subInter)-1)=subInter;
            additionalSlots=additionalSlots+numel(subInter)-1;
        else
            sub(jj)=subInter;
            
        end
    end
end

function [vectorEntries]=RemoveIdenticalEntries(vectorEntries)
    % Function which removes identical entries in a column vector
    vectorEntriesUnsort=vectorEntries;
    [vectorEntries,vectorIndex]=sort(vectorEntries);
    kk=1;
    rmvDI=[];
    for ii=2:length(vectorEntries)
        if vectorEntries(ii)==vectorEntries(ii-1)
            rmvDI(kk)=ii;
            kk=kk+1;
        end
    end
    %vectorEntries(rmvDI)=[];
    vectorIndex(rmvDI)=[];
    vectorEntries=vectorEntriesUnsort(vectorIndex);
end

%}

