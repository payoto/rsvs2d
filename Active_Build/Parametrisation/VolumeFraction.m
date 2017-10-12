%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%                - Snakes -
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Volume fraction calculation

function [volumefraction,coeffstruct,cellCentredGrid]=VolumeFraction(snaxel,snakposition,refinedGrid,volfracconnec,...
        cellCentredGrid,insideContourInfo)
    % Calculates teh volume fraction in the old cells
    
    
    insideContourInfoIndex=[refinedGrid.edge(logical(insideContourInfo)).index];
    %     [insideContourInfoIndex]=ReferenceCompArray(refinedGrid.edge,logical(insideContourInfo),'log','index');
    [cellCentredGrid]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredGrid,snakposition);
    
    for ii=1:length(cellCentredGrid)
        [cellCentredGrid(ii).areaBlock,cellCentredGrid(ii).filledvolume,cellCentredGrid(ii).coeffBlock]=...
            ExtractCellFillInformation(cellCentredGrid(ii),insideContourInfoIndex);
    end
    [volumefraction]=ExtractVolumeFractions(cellCentredGrid,volfracconnec);
    %[coeffstruct,cellCentredGrid]=VolumeFractionDifferentiated(cellCentredGrid);
    coeffstruct=[cellCentredGrid(:).coeffBlock];
    for ii=1:length(coeffstruct)
        
        coeffstruct(ii).value=CalculateVolDerivCoeff(coeffstruct(ii));
    end
end

function [volumefraction]=ExtractVolumeFractions(cellCentredGrid,volfracconnec)
    % calculates the volume fraction for each cell of the old grid
    
    volumefraction=volfracconnec.cell;
    newEdgeIndices=[cellCentredGrid(:).index];
    
    for ii=1:length(volumefraction)
        
        
        newSubs=FindObjNum([],...
            volumefraction(ii).newCellInd,newEdgeIndices);
        volumefraction(ii).splitfraction=[cellCentredGrid(newSubs).filledvolume];
        volumefraction(ii).splitvolume=[cellCentredGrid(newSubs).volume];
        volumefraction(ii).totalfraction=sum(volumefraction(ii).splitfraction);
        volumefraction(ii).totalvolume=sum(volumefraction(ii).splitvolume);
        volumefraction(ii).volumefraction=volumefraction(ii).totalfraction/...
            volumefraction(ii).totalvolume;
        if volumefraction(ii).volumefraction-1>10^-12 || volumefraction(ii).volumefraction<-10^-12
            warning('volumefraction>1')
        end
        volumefraction(ii).isSnax=false;
        jj=1;
        nSearch=length(newSubs);
        while (~volumefraction(ii).isSnax && jj<=nSearch)
            
            volumefraction(ii).isSnax=volumefraction(ii).isSnax  ...
                || ~isempty(cellCentredGrid(newSubs(jj)).snaxel);
            
            jj=jj+1;
        end
        
    end
    
    
end

function [cellCentredGrid]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredGrid,snakposition)
    % Extracts the snaxel data and matches it to the cells
    
    cellCentredGrid(1).snaxel=struct([]);
    
    %     [snakPosInd]=ReferenceCompArray(snakposition,inf,'inf','index');
    %     [edgeInd]=ReferenceCompArray(refinedGrid.edge,inf,'inf','index');
    %     [cellInd]=ReferenceCompArray(cellCentredGrid,inf,'inf','index');
    %     [vertIndex]=ReferenceCompArray(refinedGrid.vertex,inf,'inf','index');
    %     [vertCoord]=ReferenceCompArrayVertCat(refinedGrid.vertex,inf,'inf','coord');
    
    snakPosInd=[snakposition(:).index];
    edgeInd=[refinedGrid.edge(:).index];
    cellInd=[cellCentredGrid(:).index];
    vertIndex=[refinedGrid.vertex(:).index];
    vertCoord=vertcat(refinedGrid.vertex(:).coord);
    
    
    
    %     snakPosFields='coord,vector,vectorprec,vectornext,normvector,edgelength';
    %     snakPosFields=VerticalStringArray(snakPosFields,',');
    snakPosFields={'coord','vector','vectorprec','vectornext','normvector','edgelength'};
    for ii=1:length(snaxel)
        snaxEdge=snaxel(ii).edge;
        snaxEdgeSub=FindObjNum(refinedGrid.edge,snaxEdge,edgeInd);
        snaxCells=refinedGrid.edge(snaxEdgeSub).cellindex;
        snaxCells=snaxCells(snaxCells~=0);
        snaxCellsSub=FindObjNum(cellCentredGrid,snaxCells,cellInd);
        snaxPosSub=FindObjNum(snakposition,snaxel(ii).index,snakPosInd);
        
        for jj=1:length(snaxCellsSub)
            iToVert=vertCoord(find(vertIndex==snaxel(ii).tovertex),:); %#ok<FNDSB> % extract vertex coordinates
            iFromVert=vertCoord(find(vertIndex==snaxel(ii).fromvertex),:); %#ok<FNDSB>
            
            snaxelCell=snaxel(ii);
            for kk=1:length(snakPosFields(:,1))
                fieldNam=deblank(snakPosFields{kk});
                snaxelCell.(fieldNam)=0;
            end
            
            snaxelCell.coord=snakposition(snaxPosSub).coord;
            snaxelCell.vector=snakposition(snaxPosSub).vector;
            snaxelCell.vectorprec=snakposition(snaxPosSub).vectorprec;
            snaxelCell.vectornext=snakposition(snaxPosSub).vectornext;
            snaxelCell.normvector=snakposition(snaxPosSub).normvector;
            snaxelCell.edgelength=norm(iToVert-iFromVert);
            
            cellCentredGrid(snaxCellsSub(jj)).snaxel=...
                [cellCentredGrid(snaxCellsSub(jj)).snaxel,snaxelCell];
        end
    end
    
end

function [areablock,volume,coeffblock]=ExtractCellFillInformation(cellStruct,insideContourInfoIndex)
    % Extracts the contour information to calculate the area
    
    nBordBlocks=length(cellStruct.snaxel)/2;
    errFlag=false;
    if nBordBlocks>0
%         [edgeSnak]=ExtractCellSnaxelConnectedPairs(nBordBlocks,cellStruct);
%         if sum(edgeSnak==0)==2
%             warning('Snaxels are in an unconventional configuration')
%             errFlag=true;
%         else
            [areablock,coeffblock]=ExtractBorderStructure(cellStruct,[],nBordBlocks);
            [volume,areablock]=CalculateCellVolume(areablock,cellStruct.volume);
            %[coeffblock]=ExtractBorderDerivStructure(cellStruct,edgeSnak,nBordBlocks);
%         end
    end
    
    if nBordBlocks<=0 || errFlag
        coeffblock=struct([]);
        areablock.border=[0,1];
        isFullCell=prod(FindObjNum([],[cellStruct.edge(:).index],...
            insideContourInfoIndex)>0);
        if isFullCell
            volume=cellStruct.volume;
        else
            volume=0;
        end
    end
    
end

function [areablock,derivblock]=ExtractBorderStructure(cellStruct,edgeSnak,nBordBlocks)
    
    snaxInd=[cellStruct.snaxel(:).index];
    snaxEdge=[cellStruct.snaxel(:).edge];
    snaxOrder=[cellStruct.snaxel(:).orderedge];
    snaxFrom=[cellStruct.snaxel(:).fromvertex];
    snaxTo=[cellStruct.snaxel(:).tovertex];
    vertList=[cellStruct.vertex(:).index];
    vertEdgeIndex=[cellStruct.edge(:).vertexindex];
    edgeList=[cellStruct.edge(:).index];
    exploredSnax=[];
    jj=0;
    snaxSubList=1:numel(snaxInd);
    while numel(exploredSnax)<numel(snaxInd)
        %ii=1;
        ii=min(setxor(snaxSubList,exploredSnax));
        snaxStart=snaxInd(ii);
        currSnaxList=ii;
        coordList=cellStruct.snaxel(ii).coord;
        exploredSnax=[exploredSnax,ii];
        currOrder=snaxOrder(ii);
        currEdge=snaxEdge(ii);
        currTo=snaxTo(ii);
        currFrom=snaxFrom(ii);
        flag=true;
        kk=0;
        while flag
            % currently will match with itself (except if vertices set ii=0)
            isEdgeActSnax=(currEdge==snaxEdge([1:ii-1,ii+1:end])) ...
                & xor(currTo>currFrom,currOrder<snaxOrder([1:ii-1,ii+1:end]));
            
            if any(isEdgeActSnax)
                % find snaxel index which is closest
                if isfinite(currOrder)
                    snaxOrderComp=abs((snaxOrder([1:ii-1,ii+1:end])-currOrder));
                else
                    snaxOrderComp=-snaxOrder([1:ii-1,ii+1:end]);
                end
                snaxOrderComp(~isEdgeActSnax)=nan;
                [~,indNewInd]=min(snaxOrderComp);
                indNewInd=indNewInd+(indNewInd>=ii && ii>0);
                newInd=snaxInd(indNewInd);
                if FindObjNum([],indNewInd,exploredSnax)~=0
                    disp('break is at 1st Break - unexpected behaviour')
                    break
                end
                % Add Snaxel
                coordList=[coordList;cellStruct.snaxel(indNewInd).coord];
                exploredSnax=[exploredSnax,indNewInd];
                currSnaxList=[currSnaxList,indNewInd];
                % follow to next snaxel along edge
                nextSnaxSub=FindObjNum([],cellStruct.snaxel(indNewInd).connectivity,snaxInd);
                nextSnaxSub=nextSnaxSub(nextSnaxSub~=0);
%                 if numel(nextSnaxSub)==2
%                     disp('Snax 2')
%                 end
                if isempty(nextSnaxSub) || numel(nextSnaxSub)>2
                    error('There was a problem trying to follow connections for area')
                end
                ii=nextSnaxSub;
                if any(sort(exploredSnax)~=unique(exploredSnax))
                    error('Error In volume fraction structure building')
                end
                if any(FindObjNum([],ii,exploredSnax)~=0)
                    %disp('break is at 2nd Break')
                    break
                end
                exploredSnax=[exploredSnax,ii];
                currSnaxList=[currSnaxList,ii];
                coordList=[coordList;cellStruct.snaxel(ii).coord];
                currOrder=snaxOrder(ii);
                currEdge=snaxEdge(ii);
                currTo=snaxTo(ii);
                currFrom=snaxFrom(ii);
                % repeat
            else
                % find vertex which matches
                coordList=[coordList;cellStruct.vertex(FindObjNum([],currFrom,vertList)).coord];
                
                currSnaxList=[currSnaxList,0];
                % Add vertex
                currEdgeSub=FindObjNum([],currEdge,edgeList);
                matchingEdgeVert=ceil(FindObjNum([],currFrom,vertEdgeIndex)/2);
                matchingEdgeVert(matchingEdgeVert==currEdgeSub)=[];
                % find next edge
                if numel(matchingEdgeVert)~=1
                    error('Problem in the contour extraction of the volume polygon')
                end
                currEdge=cellStruct.edge(matchingEdgeVert).index;
                
                currTo=currFrom;
                currFrom=cellStruct.edge(matchingEdgeVert).vertexindex;
                currFrom(currFrom==currTo)=[];
                if currTo>currFrom
                    currOrder=Inf;
                else
                    currOrder=0;
                end
                ii=0;
                % look for snaxel on next edge
                % (currOrder=0 if initial vertex, currOrder=Inf if initial is final)
                % currTo==itself; currFrom==nextVertex
                % ii=0
                kk=kk+1;
                if kk>4
                    %warning('Something could be wrong')
                end
            end
            flag=true;
        end
        jj=jj+1;
        [areablock(jj).blockstruct]=BuildAreaBlock(coordList,cellStruct.snaxel,currSnaxList);
        [derivblock(jj).blockstruct]=BuildAreaBlockDeriv(coordList,...
            currSnaxList,cellStruct.snaxel,cellStruct.index);
        
    end
    derivblock=[derivblock(:).blockstruct];
    
end

function [areaBlock]=BuildAreaBlock(coordList,cellSnax,currSnaxList)
    
   % if ~CCWLoop(RemoveIdenticalConsecutivePoints(coordList))
   ccwtest=~CCWLoop(coordList);
   if isempty(ccwtest)
       ccwtest=false;
   end
   
   if ~CCWLoopSnax(cellSnax,currSnaxList) %|| ccwtest
        coordList=flip(coordList);
    end
    n=size(coordList,1);
    areaBlock=repmat(struct('length',[],'centre',[0 0],'normal',[0 0]),[1,n]);
    rotCW=[0 1;-1 0];
    for ii=1:n
        iip1=mod(ii,n)+1;
        areaBlock(ii).length=1;%sqrt(sum((coordList(ii,:)-coordList(iip1,:)).^2));
        areaBlock(ii).centre=(coordList(ii,:)+coordList(iip1,:))/2;
        areaBlock(ii).normal=(rotCW*(coordList(iip1,:)-coordList(ii,:))')';
    end
    
end

function [isCCW]=CCWLoopSnax(cellSnax,currSnaxList)
    
    n=length(currSnaxList);
    isCCW=[];
    for ii=n:-1:2
        iip1=mod(ii,n)+1;
        if all(currSnaxList([ii,iip1])~=0)
            
            if cellSnax(currSnaxList(ii)).snaxnext==cellSnax(currSnaxList(iip1)).index
                isCCW=true;
                break
            elseif cellSnax(currSnaxList(ii)).snaxprec==cellSnax(currSnaxList(iip1)).index
                isCCW=false;
                break
            else
                warning('should not get here.')
            end
        end
        
    end
    
    if isempty(isCCW);error('CCW direction unset');end
end

function [edgeSnak]=ExtractCellSnaxelConnectedPairs(nBordBlocks,cellStruct)
    % Extract the connected snaxels from cell information
    snaxNext=vertcat(cellStruct.snaxel(:).snaxnext);
    snaxInd=[cellStruct.snaxel(:).index];
    
    edgeSnak=zeros([nBordBlocks,2]);
    kk=0;
    for ii=1:length(snaxInd)
        nextIsInCell=sum(snaxNext(ii)==snaxInd)>0;
        if nextIsInCell
            kk=kk+1;
            edgeSnak(kk,:)=[snaxInd(ii),snaxNext(ii)];
            
        end
    end
    
    % test statements
    testEdgeSnak=numel(RemoveIdenticalEntries(edgeSnak(:)))...
        ~=numel(RemoveIdenticalEntries(snaxInd));
    if kk~=nBordBlocks || testEdgeSnak
        warning('You''re fucking up')
    end
    
end

function [volume]=CalculateGreensVolume(bordstruct)
    % THis function uses green's theorem to calculate the area inside a set
    % of snaxels
    volumeArray=zeros(size(bordstruct));
    for ii=1:length(bordstruct)
        volumeArray(ii)=dot(bordstruct(ii).centre,bordstruct(ii).normal)...
            *bordstruct(ii).length;
    end
    volume=sum(volumeArray)/2;
    
end

function [volume,areablock]=CalculateCellVolume(areablock,totalVol)
    % Calculate the volume bounded by snaxels within a cell
    
    for ii=1:length(areablock)
        areablock(ii).volume=CalculateGreensVolume(areablock(ii).blockstruct);
    end
    volume=sum(([areablock(:).volume]));
    
    if (volume-totalVol)>1e-12
        volume=totalVol+sum(([areablock(:).volume]-totalVol));
    end
    
    if volume<0
        warning('Volume<0')
%         volume=abs(volume);
%         for ii=1:length(areablock)
%             areablock(ii).blockstruct=areablock(ii).blockstruct;
%         end
    end
end

%{
function [posDelim2]=regexpREPLACE(str,delim)
    
    lD=length(delim);
    lStr=length(str);
    posDelim=zeros(size(str));
    kk=1;
    for ii=1:(lStr-(lD-1))
        
        if strcmp(delim,str(ii:(ii+lD-1)))
            posDelim(kk)=ii;
            kk=kk+1;
        end
        
    end
    posDelim2=posDelim(posDelim~=0);
    
end

function strVert=VerticalStringArray(str,delim)
    
    delimPos=regexp(str,delim);
    wS=[1,delimPos+1];
    wE=[delimPos,length(str)];
    wL=wE-wS;
    
    maxWL=max(wL);
    
    strVert=repmat(blanks(maxWL),[length(wL),1]);
    for ii=1:length(wL)
        strVert(ii,1:wL(ii))=str(wS(ii):wE(ii)-1);
    end
    
    
end

function [bordstruct,edgeSnakSub]=AddSnaxelBorders(cellStruct,edgeSnakSub)
    % Calculates snaxel-snaxel border
    % in CCW order
    
    isCCWOrder=cellStruct.snaxel(edgeSnakSub(1)).snaxnext==...
        cellStruct.snaxel(edgeSnakSub(2)).index;
    if ~isCCWOrder
        edgeSnakSub=edgeSnakSub(2:-1:1);
    end
    
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    
    normalVector=cellStruct.snaxel(edgeSnakSub(1)).vectornext;
    normNormalVec=norm(normalVector);
    if normNormalVec~=0
        normalVector=normalVector/norm(normalVector);
    end
    bordstruct=BorderStructure(norm(activeCoord(1,:)-activeCoord(2,:)),...
        mean(activeCoord),normalVector);
    [edgeSnakSub]=CalculateCWdirectionEdge(activeCoord,...
        normalVector,edgeSnakSub);
end

function [bordStruct,activeFromVertex,activeSnaxEdge]=...
        AddEdgeSnaxelBorders(cellStruct,snakSub,vertInd,edgeInd)
    % Creates teh Border structure for Edge snaxel pairs
    
    activeSnaxEdge=cellStruct.snaxel(snakSub).edge;
    activeSnaxEdgeSub=FindObjNum(cellStruct.edge,activeSnaxEdge,edgeInd);
    activeFromVertex=cellStruct.snaxel(snakSub).fromvertex;
    activeFromVertexSub=FindObjNum(cellStruct.vertex,activeFromVertex,vertInd);
    
    actCoord=zeros(size([2,length(cellStruct.snaxel(snakSub).coord)]));
    
    actCoord(1,:)=cellStruct.snaxel(snakSub).coord;
    actCoord(2,:)=cellStruct.vertex(activeFromVertexSub).coord;
    
    bordStruct=BorderStructure(norm(actCoord(1,:)-actCoord(2,:)),...
        mean(actCoord),cellStruct.edge(activeSnaxEdgeSub).normalvector);
    
    
    
end

function [bordstruct]=AddEdgeBorders(activeStartVertex,previousEdge,cellStruct,edgeSnakSub...
        ,edgeInd,vertEdgeInd)
    % Adds the borders which are related to edges
    
    bordstructTemp=struct('length',[],'centre',[0 0],'normal',[0 0]);
    bordstruct=repmat(bordstructTemp,[1 length(cellStruct.edge)]);
    
    ll=1;
    while activeStartVertex~=cellStruct.snaxel(edgeSnakSub).fromvertex
        
        previousEdgeSub=FindObjNum(cellStruct.edge,previousEdge,edgeInd);
        
        
        nextEdgeLog=logical(sum(vertEdgeInd==activeStartVertex,2));
        nextEdgeLog(previousEdgeSub)=false;
        nextEdge=cellStruct.edge(nextEdgeLog).index;
        
        bordstruct(ll)=BorderStructure(...
            cellStruct.edge(nextEdgeLog).edgelength,...
            cellStruct.edge(nextEdgeLog).edgecentre,...
            cellStruct.edge(nextEdgeLog).normalvector);
        ll=ll+1;
        
        
        nextEdgeVertices=cellStruct.edge(nextEdgeLog).vertexindex;
        nextActiveVertex=nextEdgeVertices(nextEdgeVertices~=activeStartVertex);
        previousEdge=nextEdge;
        activeStartVertex=nextActiveVertex;
    end
    bordstruct(ll:end)=[];
end

function [edgeSnak]=ExtractCellSnaxelConnectedPairsProblem(nBordBlocks,cellStruct)
    % Extract the connected snaxels from cell information
    snaxConn=vertcat(cellStruct.snaxel(:).connectivity);
    snaxInd=[cellStruct.snaxel(:).index];
    snaxIndWorking=snaxInd;
    snaxConnWorking=snaxConn;
    edgeSnak=zeros([nBordBlocks,2]);
    for ii=1:nBordBlocks
        connLog=false(size(snaxIndWorking));
        for jj=1:2
            connLog=connLog | (snaxConnWorking(1,jj)==snaxIndWorking);
        end
        connLog(1)=true;
        edgeSnak(ii,:)=snaxIndWorking(connLog)';
        snaxIndWorking(connLog)=[];
        connRmvSub=find(connLog);
        snaxConnWorking(connRmvSub,:)=[]; %#ok<FNDSB>
    end
    if ~isempty(snaxIndWorking)
        warning('You''re fucking up')
    end
    
end


function [edgeSnakSub]=CalculateCWdirectionEdge(activeCoord,normalVector,edgeSnakSub)
    
    testVector=activeCoord-(ones(2,1)*mean(activeCoord));
    baseVector=normalVector;
    [vecAngles]=ExtractAngle360(baseVector,testVector);
    [~,orderedInd]=sort(vecAngles);
    edgeSnakSub=edgeSnakSub(orderedInd);
end

function [bordstruct]=BorderStructure(bordLength,bordCentre,bordNormal)
    
    
    bordstruct.length=bordLength;
    bordstruct.centre=bordCentre;
    bordstruct.normal=bordNormal;
    
end
%}

%% Volume Derivative calculation


function [areaBlock]=BuildAreaBlockDeriv(coordList,currSnaxList,...
        cellsnax,cellInd)
    
    %if ~CCWLoop(RemoveIdenticalConsecutivePoints(coordList))
    if ~CCWLoopSnax(cellsnax,currSnaxList)
        coordList=flip(coordList);
        currSnaxList=flip(currSnaxList);
    end
    n=size(coordList,1);
    areaBlock=repmat(struct('cellindex',[],...
        'snaxelindex',[],...
        'diffgrid',[],...
        'centre',[],...
        'normal',[],...
        'ordercorner',[]),[1,sum(currSnaxList~=0)*2]);
    rotCW=[0 1;-1 0];
    kk=1;
    for ii=1:n
        iip1=mod(ii,n)+1;
        if currSnaxList(ii)>0
            areaBlock(kk).centre=(coordList(ii,:)+coordList(iip1,:))/2;
            areaBlock(kk).normal=(rotCW*(coordList(iip1,:)-coordList(ii,:))')';
            areaBlock(kk).ordercorner=1;
            areaBlock(kk).diffgrid=cellsnax(currSnaxList(ii)).vector*cellsnax(currSnaxList(ii)).edgelength;
            areaBlock(kk).snaxelindex=cellsnax(currSnaxList(ii)).index;
            areaBlock(kk).cellindex=cellInd;
            kk=kk+1;
        end
        if currSnaxList(iip1)>0
            areaBlock(kk).centre=(coordList(ii,:)+coordList(iip1,:))/2;
            areaBlock(kk).normal=(rotCW*(coordList(iip1,:)-coordList(ii,:))')';
            areaBlock(kk).ordercorner=2;
            areaBlock(kk).diffgrid=cellsnax(currSnaxList(iip1)).vector*cellsnax(currSnaxList(iip1)).edgelength;
            areaBlock(kk).snaxelindex=cellsnax(currSnaxList(iip1)).index;
            areaBlock(kk).cellindex=cellInd;
            kk=kk+1;
        end
    end
    
end

function [value]=CalculateVolDerivCoeff(coeffstruct)
    % calculates the value of the velocity derivative coefficient
    
    DgVec=coeffstruct.diffgrid;
    normalVec=coeffstruct.normal;
    positionVec=coeffstruct.centre;
    orderInd= coeffstruct.ordercorner;
    rotNeg90=[0 -1; 1 0];
    
    value=(dot(DgVec,normalVec)/2)...
        +(((-1)^(orderInd-1))*dot((rotNeg90*DgVec'),positionVec));
    
end

%{
function [areablock]=ExtractBorderDerivStructureLegacy(cellStruct,edgeSnak,nBordBlocks)
    % Extracts the informaion for the calculation of the derivative of the
    % Area
    vertInd=[cellStruct.vertex(:).index];
    edgeInd=[cellStruct.edge(:).index];
    snaxInd=[cellStruct.snaxel(:).index];
    
    edgeSnakSub=FindObjNum([],edgeSnak(:,1),snaxInd);
    edgeSnakSub(:,2)=FindObjNum([],edgeSnak(:,2),snaxInd);
    areablockTemp=struct('areablock',struct('length',[],'centre',[0 0],'normal',[0 0]));
    areablock=repmat(areablockTemp,[1 nBordBlocks]);
    for ii=nBordBlocks:-1:1
        % Calculate snaxel to snaxel border
        [bordstruct1,edgeSnakSub(ii,:)]=AddSnaxelDeriv(cellStruct,edgeSnakSub(ii,:));
        % First snax-edge border
        [bordstruct2,~,~]=...
            AddEdgeSnaxelDeriv(cellStruct,edgeSnakSub(ii,1),vertInd,edgeInd,2);
        % First snax-edge border
        [bordstruct3,~,~]=...
            AddEdgeSnaxelDeriv(cellStruct,edgeSnakSub(ii,2),vertInd,edgeInd,1);
        
        areablock(ii).blockstruct=[bordstruct1,bordstruct2,bordstruct3];
    end
    areablock=[areablock(:).blockstruct];
    
end

function [areablock]=ExtractBorderDerivStructure(cellStruct,edgeSnak,nBordBlocks)
    
    snaxInd=[cellStruct.snaxel(:).index];
    snaxEdge=[cellStruct.snaxel(:).edge];
    snaxOrder=[cellStruct.snaxel(:).orderedge];
    snaxFrom=[cellStruct.snaxel(:).fromvertex];
    snaxTo=[cellStruct.snaxel(:).tovertex];
    vertList=[cellStruct.vertex(:).index];
    vertEdgeIndex=[cellStruct.edge(:).vertexindex];
    edgeList=[cellStruct.edge(:).index];
    exploredSnax=[];
    jj=0;
    snaxSubList=1:numel(snaxInd);
    
    
    while numel(exploredSnax)<numel(snaxInd)
        %ii=1;
        
        ii=min(setxor(snaxSubList,exploredSnax));
        snaxStart=snaxInd(ii);
        currSnaxList=ii;
        coordList=cellStruct.snaxel(ii).coord;
        exploredSnax=[exploredSnax,ii];
        currOrder=snaxOrder(ii);
        currEdge=snaxEdge(ii);
        currTo=snaxTo(ii);
        currFrom=snaxFrom(ii);
        flag=true;
        kk=0;
        while flag
            % currently will match with itself (except if vertices set ii=0)
            isEdgeActSnax=(currEdge==snaxEdge([1:ii-1,ii+1:end])) ...
                & xor(currTo>currFrom,currOrder<snaxOrder([1:ii-1,ii+1:end]));
            
            if any(isEdgeActSnax)
                % find snaxel index which is closest
                if isfinite(currOrder)
                    snaxOrderComp=abs((snaxOrder([1:ii-1,ii+1:end])-currOrder));
                else
                    snaxOrderComp=snaxOrder([1:ii-1,ii+1:end]);
                end
                snaxOrderComp(~isEdgeActSnax)=nan;
                [~,indNewInd]=min(snaxOrderComp);
                indNewInd=indNewInd+(indNewInd>=ii && ii>0);
                newInd=snaxInd(indNewInd);
                if FindObjNum([],indNewInd,exploredSnax)~=0
                    disp('break is at 1st Break - unexpected behaviour')
                    break
                end
                % Add Snaxel
                coordList=[coordList;cellStruct.snaxel(indNewInd).coord];
                exploredSnax=[exploredSnax,indNewInd];
                currSnaxList=[currSnaxList,indNewInd];
                % follow to next snaxel along edge
                nextSnaxSub=FindObjNum([],cellStruct.snaxel(indNewInd).connectivity,snaxInd);
                nextSnaxSub=nextSnaxSub(nextSnaxSub~=0);
                if isempty(nextSnaxSub) || numel(nextSnaxSub)>2
                    error('There was a problem trying to follow connections for area')
                end
                ii=nextSnaxSub;
                if any(sort(exploredSnax)~=unique(exploredSnax))
                    error('Error In volume fraction structure building')
                end
                if FindObjNum([],ii,exploredSnax)~=0
                    %disp('break is at 2nd Break')
                    break
                end
                exploredSnax=[exploredSnax,ii];
                currSnaxList=[currSnaxList,ii];
                coordList=[coordList;cellStruct.snaxel(ii).coord];
                currOrder=snaxOrder(ii);
                currEdge=snaxEdge(ii);
                currTo=snaxTo(ii);
                currFrom=snaxFrom(ii);
                % repeat
            else
                % find vertex which matches
                coordList=[coordList;cellStruct.vertex(FindObjNum([],currFrom,vertList)).coord];
                
                currSnaxList=[currSnaxList,0];
                % Add vertex
                currEdgeSub=FindObjNum([],currEdge,edgeList);
                matchingEdgeVert=ceil(FindObjNum([],currFrom,vertEdgeIndex)/2);
                matchingEdgeVert(matchingEdgeVert==currEdgeSub)=[];
                % find next edge
                if numel(matchingEdgeVert)~=1
                    error('Problem in the contour extraction of the volume polygon')
                end
                currEdge=cellStruct.edge(matchingEdgeVert).index;
                
                currTo=currFrom;
                currFrom=cellStruct.edge(matchingEdgeVert).vertexindex;
                currFrom(currFrom==currTo)=[];
                if currTo>currFrom
                    currOrder=Inf;
                else
                    currOrder=0;
                end
                ii=0;
                % look for snaxel on next edge
                % (currOrder=0 if initial vertex, currOrder=Inf if initial is final)
                % currTo==itself; currFrom==nextVertex
                % ii=0
                kk=kk+1;
                if kk>4
                    error('Something is probably wrong')
                end
            end
            flag=true;
        end
        jj=jj+1;
        [areablock(jj).blockstruct]=BuildAreaBlockDeriv(coordList,...
            currSnaxList,cellStruct.snaxel,cellStruct.index);
        
    end
    areablock=[areablock(:).blockstruct];
end

function [coeffstruct,cellCentredGrid]=VolumeFractionDifferentiated(cellCentredGrid)
    % Calculates the volume coefficients for the derivative of the volume fraction
    % In order to calculate the velocities
    
    for ii=1:length(cellCentredGrid)
        [cellCentredGrid(ii).coeffBlock]=...
            ExtractCellDerivInformation(cellCentredGrid(ii));
    end
    
    coeffstruct=[cellCentredGrid(:).coeffBlock];
    for ii=1:length(coeffstruct)
        
        coeffstruct(ii).value=CalculateVolDerivCoeff(coeffstruct(ii));
    end
    
end

function [areablock]=ExtractBorderStructureLegacy(cellStruct,edgeSnak,nBordBlocks)
    % Extracts the border information and sets it out in the right
    % structure form
    vertInd=[cellStruct.vertex(:).index];
    edgeInd=[cellStruct.edge(:).index];
    vertEdgeInd=vertcat(cellStruct.edge(:).vertexindex);
    
    snaxInd=[cellStruct.snaxel(:).index];
    edgeSnakSub=FindObjNum([],edgeSnak(:,1),snaxInd);
    edgeSnakSub(:,2)=FindObjNum([],edgeSnak(:,2),snaxInd);
    areablockTemp=struct('areablock',struct('length',[],'centre',[0 0],'normal',[0 0]));
    areablock=repmat(areablockTemp,[1 nBordBlocks]);
    for ii=nBordBlocks:-1:1
        
        % Calculate snaxel to snaxel border
        [bordstruct1,edgeSnakSub(ii,:)]=AddSnaxelBorders(cellStruct,edgeSnakSub(ii,:));
        % First snax-edge bord
        [bordstruct2,activeFromVertex,activeSnaxEdge]=...
            AddEdgeSnaxelBorders(cellStruct,edgeSnakSub(ii,1),vertInd,edgeInd);
        % Edge Borders
        activeStartVertex=activeFromVertex;
        previousEdge=activeSnaxEdge;
        bordstructEdges=AddEdgeBorders(activeStartVertex,previousEdge,...
            cellStruct,edgeSnakSub(ii,2),edgeInd,vertEdgeInd);
        
        % First snax-edge bord
        [bordstruct3,~,~]=...
            AddEdgeSnaxelBorders(cellStruct,edgeSnakSub(ii,2),vertInd,edgeInd);
        areablock(ii).blockstruct=[bordstruct1,bordstruct2,bordstructEdges,bordstruct3];
    end
    
end

function [coeffblock]=ExtractCellDerivInformation(cellStruct)
    % Extracts the contour information to calculate the area
    
    nBordBlocks=length(cellStruct.snaxel)/2;
    errFlag=false;
    if nBordBlocks>0
        [edgeSnak]=ExtractCellSnaxelConnectedPairs(nBordBlocks,cellStruct);
        if sum(edgeSnak==0)==2
            warning('Snaxels are in an unconventional configuration')
            errFlag=true;
        end
        if edgeSnak~=0
            [coeffblock]=ExtractBorderDerivStructure(cellStruct,edgeSnak,nBordBlocks);
        end
    end
    
    if nBordBlocks<=0 || errFlag
        coeffblock=struct([]);
    end
    
end

function [derivstruct,edgeSnakSub]=AddSnaxelDeriv(cellStruct,edgeSnakSub)
    % Calculates snaxel-snaxel border
    
    isCCWOrder=cellStruct.snaxel(edgeSnakSub(1)).snaxnext==cellStruct.snaxel(edgeSnakSub(2)).index;
    if ~isCCWOrder
        edgeSnakSub=edgeSnakSub(2:-1:1);
    end
    
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    tanVector=activeCoord(2,:)-activeCoord(1,:);
    normalVector=cellStruct.snaxel(edgeSnakSub(1)).vectornext;
    
    activeCoord=vertcat(cellStruct.snaxel(edgeSnakSub).coord);
    normVecLength=normalVector*norm(tanVector);
    
    derivstructTemp=struct('cellindex',[],...
        'snaxelindex',[],...
        'diffgrid',[],...
        'centre',[],...
        'normal',[],...
        'ordercorner',[]);
    derivstruct=repmat(derivstructTemp,[1 2]);
    
    for ii=2:-1:1
        jj=ii;%abs(ii-3);
        derivstruct(ii)=...
            BorderDerivCoeffStructure(cellStruct.index,cellStruct.snaxel(edgeSnakSub(ii)).index,...
            cellStruct.snaxel(edgeSnakSub(ii)).vector*cellStruct.snaxel(edgeSnakSub(ii)).edgelength...
            ,mean(activeCoord),normVecLength,jj);
    end
    
    
end

function [derivStruct,activeFromVertex,activeSnaxEdge]=...
        AddEdgeSnaxelDeriv(cellStruct,snakSub,vertInd,edgeInd,snaxelBorderOrder)
    % Creates the Border structure for Edge snaxel pairs
    
    activeSnaxEdge=cellStruct.snaxel(snakSub).edge;
    activeSnaxEdgeSub=FindObjNum(cellStruct.edge,activeSnaxEdge,edgeInd);
    activeFromVertex=cellStruct.snaxel(snakSub).fromvertex;
    activeFromVertexSub=FindObjNum(cellStruct.vertex,activeFromVertex,vertInd);
    
    actCoord=zeros([2,length(cellStruct.snaxel(snakSub).coord)]);
    
    actCoord(1,:)=cellStruct.snaxel(snakSub).coord;
    actCoord(2,:)=cellStruct.vertex(activeFromVertexSub).coord;
    
    normVecLength=cellStruct.edge(activeSnaxEdgeSub).normalvector...
        *norm(actCoord(1,:)-actCoord(2,:));
    
    
    derivStruct=...
        BorderDerivCoeffStructure(cellStruct.index,cellStruct.snaxel(snakSub).index,...
        cellStruct.snaxel(snakSub).vector*cellStruct.snaxel(snakSub).edgelength,...
        mean(actCoord),normVecLength,snaxelBorderOrder);
    
end

function [derivstruct]=BorderDerivCoeffStructure(cellIndex,snaxelIndex,Dg,bordCentre,bordNormal,ordercorner)
    
    derivstruct.cellindex=cellIndex;
    derivstruct.snaxelindex=snaxelIndex;
    derivstruct.diffgrid=Dg;
    derivstruct.centre=bordCentre;
    derivstruct.normal=bordNormal;
    derivstruct.ordercorner=ordercorner; % Defines wether is negative or positive in coeff cal
    
end
%}
%% Compilation Utilities

function [array]=ReferenceCompArray(obj,ref,typeRef,fields)
    
    if strcmp(typeRef,'inf')
        array=zeros([1 length(obj)]);
        for ii=1:length(array)
            array(ii)=obj(ii).(fields);
        end
    elseif strcmp(typeRef,'log')
        array=zeros([1 sum(ref)]);
        kk=1;
        for ii=1:length(ref)
            if ref(ii)
                array(kk)=obj(ii).(fields);
                kk=kk+1;
            end
        end
    elseif strcmp(typeRef,'ind')
        array=zeros([1 length(ref)]);
        for ii=ref
            
            array(kk)=obj(ii).(fields);
            kk=kk+1;
            
        end
    else
        error('Not a valid type');
    end
    
    
end

function [array]=ReferenceCompArrayVertCat(obj,ref,typeRef,fields)
    
    nCat=length(obj(1).(fields));
    
    if strcmp(typeRef,'inf')
        array=zeros([length(obj),nCat]);
        for ii=1:length(array)
            for jj=1:nCat
                array(ii,jj)=obj(ii).(fields)(jj);
            end
        end
    elseif strcmp(typeRef,'log')
        array=zeros([sum(ref) nCat]);
        kk=1;
        for ii=1:length(array)
            if ref(ii)
                for jj=1:nCat
                    array(kk,jj)=obj(ii).(fields)(jj);
                end
                kk=kk+1;
            end
        end
    elseif strcmp(typeRef,'ind')
        array=zeros([length(ref),nCat]);
        for ii=ref
            for jj=1:nCat
                array(kk,jj)=obj(ii).(fields)(jj);
            end
            kk=kk+1;
            
        end
    else
        error('Not a valid type');
    end
    
    
end