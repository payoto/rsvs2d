%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%          - Grid Refinement code -
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [gridrefined2,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement(gridreshape,param)
    % bridge function for the refinement function
    
    %unpacking input parameters
    varExtract={'refineGrid','boundstr','typeRefine','execTest','cellRefineShape','cellGeometry'};
    [refineGrid,boundstr,typeRefine,execTest,cellRefineShape,cellGeometry]=...
        ExtractVariables(varExtract,param);
    nRefine=refineGrid;
    if execTest
        ExecuteTestCode();
    end
    
    if strcmp(cellGeometry,'triangle') && ...
            ~strcmp(cellRefineShape,'fullytriangle')
        cellRefineShape='fullytriangle';
        warning('Incompatible cellGeometry and cellRefineShape - defaulting to fullytriangle')
    elseif ~strcmp(cellGeometry,'triangle') && ...
            strcmp(cellRefineShape,'fullytriangle')
        cellRefineShape='triangle';
        warning('Incompatible cellGeometry and cellRefineShape - defaulting to triangle')
    end
    
    
    switch cellRefineShape
        case 'square'
            [gridrefined,connectstructinfo]=RefineGrid(gridreshape,nRefine,typeRefine);
        case 'triangle'
            [gridrefined,connectstructinfo]=RefineGrid(gridreshape,nRefine(unique([1,1:end-1])),typeRefine);
            [gridrefined,connectstructinfo]=MakeTriangularGrid(gridrefined,connectstructinfo,nRefine(end));
        case 'fullytriangle'
            [gridrefined,connectstructinfo]=RefineTriangular(gridreshape,nRefine(1));
    end
    [gridrefined2]=EdgePropertiesReshape(gridrefined);
    [loop]=GenerateSnakStartLoop(gridrefined2,boundstr);
    
    [unstructuredrefined]=ModifReshape(gridrefined2);
end


function [loop]=GenerateSnakStartLoop(gridrefined2,boundstr)
    
    isEdge=[gridrefined2.edge(:).(boundstr{1})];
    cond=boundstr{3};
    [loop]=OrderSurfaceVertexReshape(gridrefined2,isEdge,cond);
    
    [loop]=EdgeInCondForVertex(loop,gridrefined2,cond);
end

function [gridrefined,connecstructinfo]=RefineGrid(gridreshape,nRefine,typeRefine)
    % Function containing the grid refinement 
    
    [gridreshape.cell(:).coord]=deal([]);
    [gridreshape.cell(:).vertexindex]=deal([]);
    % Pre refinement processing
    [cellRefine,cellRefinePos]=IdentifyRefineCell(gridreshape,typeRefine);
    if numel(nRefine)==1
        [gridrefined,connecstructinfo]=GridRefine_Wrapper(gridreshape,[cellRefine;cellRefinePos]',[nRefine nRefine]);
    else
        [gridrefined,connecstructinfo]=GridRefine_Wrapper(gridreshape,[cellRefine;cellRefinePos]', nRefine(1:2));
    end
    
    % Shit codfe written under the influence
    connecstructinfoCell=connecstructinfo.cell;
    kk=1;
    for ii=1:numel(gridreshape.cell)
        connecstructinfo.cell(ii).old=gridreshape.cell(ii).index;
        if gridreshape.cell(ii).index==connecstructinfoCell(kk).old;
            connecstructinfo.cell(ii).new=connecstructinfoCell(kk).new;
            kk=kk+1;
            kk=min(kk,numel(connecstructinfoCell));
        else
        connecstructinfo.cell(ii).new=gridreshape.cell(ii).index;
        end
    end
    
    for ii=1:numel(gridrefined.cell)
        gridrefined.cell(ii).refineVec(gridrefined.cell(ii).refineVec==0)=1;
    end
    
    oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
        {connecstructinfo.cell(:).new},...
        {connecstructinfo.cell(:).old},'UniformOutput',false));
    newSubs=FindObjNum([],[connecstructinfo.cell(:).new],[gridrefined.cell(:).index]);
    oldSubs=FindObjNum([],oldIndsNewOrd,[gridreshape.cell(:).index]);
    
    [gridrefined.cell(newSubs).isactive]=deal(gridreshape.cell(oldSubs).isactive);
    
%     CheckResultsRefineGrid(gridreshape)
%     CheckResultsRefineGrid(gridedgerefined)
%     CheckResultsRefineGrid(gridrefined)
    
end

%% Refinement preparation operations

function [cellRefine,cellRefinePos]=IdentifyRefineCell(gridreshape,typeRefine)
    % Identifies which cell need to be refined based on the type of
    % refinement call required
    % Returns the indices of the cells that must be refined
   
    switch typeRefine
        case 'all'
            cellRefine=[gridreshape.cell(:).index];
            cellRefinePos=1:length(gridreshape.cell);
        case 'grey'
            cellRefineLog=[gridreshape.cell(:).fill]~=1 & [gridreshape.cell(:).fill]~=0;
            cellRefine=[gridreshape.cell(cellRefineLog).index];
            cellRefinePos=find(cellRefineLog);
        case 'actgrey'
            cellRefineLog=([gridreshape.cell(:).fill]~=1 & [gridreshape.cell(:).fill]~=0)...
                | logical([gridreshape.cell(:).isactive]);
            cellRefine=[gridreshape.cell(cellRefineLog).index];
            cellRefinePos=find(cellRefineLog);
            
        case 'active'
            cellRefineLog=logical([gridreshape.cell(:).isactive]);
            cellRefine=[gridreshape.cell(cellRefineLog).index];
            cellRefinePos=find(cellRefineLog);
        case 'special'
            cellRefine=[20 21];
            cellRefinePos=[20,21];
        case 'automatic'
            cellRefineLog=logical([gridreshape.cell(:).isrefine]);
            cellRefine=[gridreshape.cell(cellRefineLog).index];
            cellRefinePos=find(cellRefineLog);
        otherwise 
            error('Unsupported refinement type')
    end

end

function [edgeRefine]=IdentifyRefineEdges(gridreshape,cellRefine)
   % Finds the edges which need to be cut to match the cells that need to
   % be refined
   
   edgeCellIndex=vertcat(gridreshape.edge(:).cellindex);
   edgeIndex=vertcat(gridreshape.edge(:).index);
   
   edgeRefineLog=zeros(size(edgeCellIndex));
   
   for ii=1:length(cellRefine)
       edgeRefineLog=edgeRefineLog+(edgeCellIndex==cellRefine(ii));
   end
   edgeRefine=edgeIndex(sum(edgeRefineLog,2)>0);
end

function [ccwinfoCell]=PrepareCellRefinement(gridreshape,cellIndex)
    % Prepares the information required for the refinement of a given cell
    
    cellIndexData=vertcat(gridreshape.edge(:).cellindex);
    correspEdgesSub=find(sum((cellIndex==cellIndexData),2));
    correspEdgesIndex=[gridreshape.edge(correspEdgesSub).index];
    correspEdgesVertex=vertcat(gridreshape.edge(correspEdgesSub).vertexindex);
    
    [cellOrderedVertex,cellOrderedEdges]=...
        OrderBlockEdges(correspEdgesVertex,...
        ones([1,length(correspEdgesVertex)])*cellIndex);
    if length(cellOrderedVertex)~=1
        error('Incorrect cell edges')
    end
    ccwinfoCell.vertex=cellOrderedVertex{1}(:,1);
    ccwinfoCell.edge=[gridreshape.edge(correspEdgesSub(cellOrderedEdges{1})).index];
    
    verticesSub=FindObjNum(gridreshape.vertex,ccwinfoCell.vertex);
    coord=vertcat(gridreshape.vertex(verticesSub).coord);
    loopStartInd=LeftMostCorner(coord);
    ccwinfoCell.vertex=ccwinfoCell.vertex([loopStartInd:end,1:loopStartInd-1]);
    ccwinfoCell.edge=ccwinfoCell.edge([loopStartInd:end,1:loopStartInd-1]);
    
    
    [isCCW]=CCWLoop(coord);
    if ~isCCW
        ccwinfoCell.vertex=ccwinfoCell.vertex([1,end:-1:2]);
        ccwinfoCell.edge=ccwinfoCell.edge(end:-1:1);
    end
    ccwinfoCell.cell=[];
end

%% Top level Refinement Operations

function [template,ccwinfo]=CreateCellRefinementTemplate(nRefine)
    % Creates the cell refinement template
    [vertextemplate,ccwinfo.vertex]=TemplateVertices(nRefine);
    [celltemplate]=TemplateCell(nRefine);
    [edgetemplate,ccwinfo.edge]=TemplateEdge(nRefine);
    ccwinfo.cell=[];
    
    template.vertex=vertextemplate;
    template.cell=celltemplate;
    template.edge=edgetemplate;
    %CheckResults(template);
end

function [gridreshape,connecstructinfo]=...
        SplitSpecifiedEdges(gridreshape,edgeIndices,nRefine)
    % Splits all the edges in the grid corresponding to a cell that will be
    % refined
    
    connecstructinfo.edge=[];
    
    nIndices=length(edgeIndices);
    
    for ii=1:nIndices
        
        [newedges,newvertices]=SplitEdge(edgeIndices(ii),nRefine,gridreshape);
        gridreshape.edge=[gridreshape.edge;newedges'];
        gridreshape.vertex=[gridreshape.vertex;newvertices'];
        connecstructinfo.edge(ii).old=edgeIndices(ii);
        connecstructinfo.edge(ii).new=[newedges(:).index];
        
    end
    
    edgeSubs=FindObjNum(gridreshape.edge,edgeIndices);
    gridreshape.edge(edgeSubs)=[];
    
end

function [gridreshape,connecstructinfo]=...
        SplitSpecifiedCells(gridreshape,cellIndices,template,ccwtempinfo)
    % Splits all the edges in the grid corresponding to a cell that will be
    % refined
    
    
    nIndices=length(cellIndices);
    
    for ii=nIndices:-1:1
            
        [ccwcellinfo]=PrepareCellRefinement(gridreshape,cellIndices(ii));
        [gridreshape,connecstructinfo(ii)]=SplitCell(template,ccwtempinfo,ccwcellinfo,gridreshape,cellIndices(ii));
        %CheckResults(gridreshape)
    end
    if ~exist('connecstructinfo','var');connecstructinfo=[];end
end

function [template,ccwinfo]=TriangularTemplate(nRefine)
    
    [template1,ccwinfo]=CreateCellRefinementTemplate(1);
    [template2,ccwinfo]=CreateCellRefinementTemplate(nRefine);
    connec.cell.oldCellInd=[1];
    connec.cell.newCellInd=[1:numel(template2.cell)];
    [template]=CoarsenGrid(template2,template1,connec);
    ii=max([template.vertex.index])+1;
    iiVert=numel(template.vertex)+1;
    template.vertex(iiVert).index=ii;
    template.vertex(iiVert).coord=[0.5 0.5];
    nEdge=numel(template.edge);
    
    edgeInd=[template.edge.index];
    for kk=1:numel(template.edge)
        template.edge(kk).index=kk;
    end
    ccwinfo.edge=FindObjNum([],ccwinfo.edge,edgeInd);
    ll=nEdge;
    edgeInd=[template.edge.index];
    for jj=1:nEdge
        kk=FindObjNum([],ccwinfo.edge(jj),edgeInd);
        cInd=template.edge(kk).cellindex;
        template.edge(kk).cellindex(cInd~=0)=jj;
        template.edge(end+1)=template.edge(1);
        ll=ll+1;
        template.edge(end).index=ll;
        template.edge(end).vertexindex=[ccwinfo.vertex(jj),ii];
        template.edge(end).cellindex=[jj,mod(jj-2,nEdge)+1];
        template.edge(end).orientation=0.5;
        template.cell(jj)=template.cell(1);
        template.cell(jj).index=jj;
    end
    
end

function [gridreshape,origconnec]=MakeTriangularGrid(gridreshape,origconnec,nRefine)
    
    [template,ccwinfo]=TriangularTemplate(nRefine);
    %[template,ccwinfo]=CreateCellRefinementTemplate(nRefine);
    
    edgeIndices=[gridreshape.edge.index];
    cellIndices=[gridreshape.cell.index];
    [gridreshape,connecstructinfo]=...
         SplitSpecifiedEdges(gridreshape,edgeIndices,nRefine);
    [gridreshape,connecstructinfo]=...
        SplitSpecifiedCells(gridreshape,cellIndices,template,ccwinfo);
    
    newOld=[connecstructinfo.old];
    for ii=1:numel(origconnec.cell)
        origconnec.cell(ii).new=[connecstructinfo(FindObjNum([],...
            origconnec.cell(ii).new,newOld)).new];
    end
    
end
%% Grid Operations

function [newedges,newvertices]=SplitEdge(edgeIndex,nRefine,gridreshape)
    % splits an edge into nRefine edges adding corresponding vertices
    newedges=AddEdgeStructure([],[],[],[]);
    newvertices=AddVertexStructure([],[]);
    edgeIndices=[gridreshape.edge(:).index];
    edgeSub=FindObjNum(gridreshape.edge,edgeIndex,edgeIndices);
    maxEdgeIndex=max(edgeIndices);
    baseedge=gridreshape.edge(edgeSub);
    
    vertexIndices=[gridreshape.vertex(:).index];
    maxVertexIndex=max(vertexIndices);
    vertexStream=[baseedge.vertexindex(1),...
        ((maxVertexIndex+1):(maxVertexIndex+nRefine-1)),...
        baseedge.vertexindex(2)];
    origVertSub=FindObjNum(gridreshape.vertex,baseedge.vertexindex,vertexIndices);
    baseCoord=vertcat(gridreshape.vertex(origVertSub).coord);
    
    for ii=nRefine:-1:1
        [newedges(ii)]=AddEdgeStructure(maxEdgeIndex+ii,baseedge.cellindex,...
            vertexStream(ii:ii+1),[]);
    end
    for ii=nRefine-1:-1:1
        coord=baseCoord(1,:)+(baseCoord(2,:)-baseCoord(1,:))*(ii/nRefine);
        [newvertices(ii)]=AddVertexStructure(vertexStream(ii+1),coord);
    end
    %vertices=[gridreshape.vertex(origVertSub),newvertices];
end

function [gridreshape,connectcellinfo]=...
        SplitCell(template,ccwtempinfo,ccwcellinfo,gridreshape,cellInd)
    % Splits a grid cell with all the corresponding operations
    
    [edgeList,vertList,cellList]=...
        CreateCellRefinementLists(template,ccwtempinfo,ccwcellinfo,gridreshape);
    [adaptedtemplate]=AdaptTemplate(template,edgeList,vertList,cellList,gridreshape,cellInd);
    
    cellDelSub=FindObjNum(gridreshape.cell,cellInd);
    edgeDelSub=FindObjNum(gridreshape.edge,edgeList);
    edgeDelSub(edgeDelSub==0)=[];
    vertDelSub=FindObjNum(gridreshape.vertex,vertList);
    vertDelSub(vertDelSub==0)=[];
    
    
    gridreshape.cell(cellDelSub)=[];
    gridreshape.edge(edgeDelSub)=[];
    gridreshape.vertex(vertDelSub)=[];
    
    gridreshape.cell=[gridreshape.cell;adaptedtemplate.cell'];
    gridreshape.vertex=[gridreshape.vertex;adaptedtemplate.vertex'];
    gridreshape.edge=[gridreshape.edge;adaptedtemplate.edge'];
    
    connectcellinfo.old=cellInd;
    connectcellinfo.new=cellList;
    
end

function [adaptedtemplate]=AdaptTemplate(template,edgeList,vertList,cellList,gridreshape,cellInd)
    % Adapts the template to match the specific indices of the case being
    % considered
   
    cornSub=FindObjNum([],vertList,[gridreshape.vertex.index]);
    cornSub=cornSub(cornSub~=0);
    baseCoord=vertcat(gridreshape.vertex(cornSub).coord);
    baseCoord=[min(baseCoord);max(baseCoord)];
    [adaptedtemplate.vertex]=AdaptTemplateVertex(template,vertList,baseCoord);
    
    cellSub=FindObjNum(gridreshape.cell,cellInd);
    fillBase=gridreshape.cell(cellSub).fill;
    [adaptedtemplate.cell]=AdaptTemplateCell(template,cellList,fillBase);
    
    [adaptedtemplate.edge]=AdaptTemplateEdge(template,vertList,edgeList,cellList,gridreshape,cellInd);
end

function [adaptedvertex]=AdaptTemplateVertex(template,vertList,baseCoord)
    % Modify Vertex Indices & Coord
    adaptedvertex=template.vertex;
    base=baseCoord(1,:);
    slope=baseCoord(2,:)-baseCoord(1,:);
    for ii=1:length(adaptedvertex)
        adaptedvertex(ii).coord=((adaptedvertex(ii).coord).*slope)+base;
        adaptedvertex(ii).index=vertList(ii);
    end
    
end

function [adaptededge]=AdaptTemplateEdge(template,vertList,edgeList,cellList,gridreshape,cellInd)

    % Modify Edges
    % index
    % cellindex
    % vertexindex
    adaptededge=template.edge;
    cellListExt=[0,cellList];
    for ii=1:length(adaptededge)
        adaptededge(ii).index=edgeList(ii);
        adaptededge(ii).cellindex=cellListExt(adaptededge(ii).cellindex+1);
        adaptededge(ii).vertexindex=vertList(adaptededge(ii).vertexindex);
    end
    
    edgeListSub=FindObjNum(gridreshape.edge,edgeList);
    edgeModif=find(edgeListSub~=0);
    
    cellIndexEdgeModif=vertcat(gridreshape.edge(edgeListSub(edgeModif)).cellindex);
    cellIndexEdgeModif(cellIndexEdgeModif==cellInd)=0;
    cellIndexEdgeKeep=sum(cellIndexEdgeModif,2);
    
    for ii=1:length(edgeModif)
        
        adaptededge(edgeModif(ii)).cellindex...
            (adaptededge(edgeModif(ii)).cellindex==0)=cellIndexEdgeKeep(ii);
    end
    
    
    
end

function [adaptedcell]=AdaptTemplateCell(template,cellList,fillBase)

    adaptedcell=template.cell;
    for ii=1:length(adaptedcell)
        adaptedcell(ii).fill=fillBase;
        adaptedcell(ii).index=cellList(ii);
    end
end

%% Simple Cell refinement steps

function [edgeList,vertList,cellList]=...
        CreateCellRefinementLists(template,ccwtempinfo,ccwcellinfo,gridreshape)
    % Matches the index of existing elements
    
    [edgeList]=ExtractObjectList...
        (template,ccwtempinfo,ccwcellinfo,gridreshape,'edge');
    [vertList]=ExtractObjectList...
        (template,ccwtempinfo,ccwcellinfo,gridreshape,'vertex');
    [cellList]=ExtractObjectList...
        (template,ccwtempinfo,ccwcellinfo,gridreshape,'cell');
    
end

function [objList]=ExtractObjectList...
        (template,ccwtempinfo,ccwcellinfo,gridreshape,objstr)
    % Extracts Object List matching template to existing geometries
    
    templateIndex=[template.(objstr)(:).index];
    nInd=length(templateIndex);
    minNewIndex=max([gridreshape.(objstr)(:).index])+1;
    ccwTempInfo=ccwtempinfo.(objstr);
    ccwCellInfo=ccwcellinfo.(objstr);
    nExist=length(ccwCellInfo);
    nNewInd=nInd-nExist;
    newIndList=minNewIndex:minNewIndex+nNewInd-1;
    objList=zeros([1,nInd]);
    objList(ccwTempInfo)=ccwCellInfo;
    objList(objList==0)=newIndList;
    
end

%% Simple template creation steps

function [vertextemplate,ccwBorderVertices]=TemplateVertices(nRefine)
    % Creates the vertex template and returns the border
    
    for ii=1:(nRefine+1)
        for jj=1:(nRefine+1)
            linInd=ii+(jj-1)*(nRefine+1);
            vertextemplate(linInd)=...
                AddVertexStructure(linInd,(([ii,jj]-1)/(nRefine)));
        end 
    end 
    ccwBorderVertices=[1:(nRefine+1-1),...
        (nRefine+1):(nRefine+1):((nRefine+1)^2-1),...
        ((nRefine+1)^2):-1:((nRefine+1)^2-nRefine+1)...
        ((nRefine+1)^2-nRefine):-(nRefine+1):(1+1)];
end

function [celltemplate]=TemplateCell(nRefine)
    % Creates the vertex template and returns the border
    
    for ii=1:(nRefine)
        for jj=1:(nRefine)
            linInd=ii+(jj-1)*(nRefine);
            [celltemplate(linInd)]=AddCellStructure(linInd,1);
            
        end 
    end 
end

function [edgetemplate,ccwBorderEdges]=TemplateEdge(nRefine)
    % Creates the edge template and returns the border
    
    [edgetemplate(1:(nRefine*(nRefine+1)))]=TemplateEdgeHorizontal(nRefine);
    [edgetemplate((1+(nRefine*(nRefine+1))):(2*(nRefine*(nRefine+1))))]=...
        TemplateEdgeVertical(nRefine);
    [ccwBorderEdges]=TemplateBorderEdges(nRefine);
    
end

function [edgetemplate]=TemplateEdgeHorizontal(nRefine)
    % Creates the template for horizontal edges
    
    for ii=1:nRefine
        for jj=1:nRefine+1
           linInd=ii+(jj-1)*(nRefine);
           indEdge=linInd;
           
           cellRight=ii+(jj-1)*(nRefine);
           cellLeft=ii+(jj-2)*(nRefine);
           cellIndex=[cellRight,cellLeft];
           cellIndex(cellIndex<1)=0;
           cellIndex(cellIndex>nRefine^2)=0;
           
           vertStart=ii+(jj-1)*(nRefine+1);
           vertEnd=(ii+1)+(jj-1)*(nRefine+1);
           vertexIndex=[vertStart,vertEnd];
           
           edgetemplate(linInd)=AddEdgeStructure(indEdge,cellIndex,vertexIndex,0);
        end
    end

    
end

function [edgetemplate]=TemplateEdgeVertical(nRefine)
    % Creates the template for vertical edges
    
    for ii=1:nRefine
        for jj=1:nRefine+1
           linInd=ii+(jj-1)*(nRefine);
           indEdge=linInd+(nRefine*(nRefine+1));
           
           cellRight=jj+(ii-1)*(nRefine);
           if jj>nRefine
               cellRight=0;
           end
           cellLeft=jj-1+(ii-1)*(nRefine);
           if jj-1<1
               cellLeft=0;
           end
           cellIndex=[cellRight,cellLeft];
           cellIndex(cellIndex<1)=0;
           cellIndex(cellIndex>nRefine^2)=0;
           
           vertStart=jj+(ii-1)*(nRefine+1);
           vertEnd=(jj)+(ii-1+1)*(nRefine+1);
           vertexIndex=[vertStart,vertEnd];
           
           edgetemplate(linInd)=AddEdgeStructure(indEdge,cellIndex,vertexIndex,1);
        end
    end
    
end

function [ccwBorderEdges]=TemplateBorderEdges(nRefine)
    % Creates the list of border edges of the template
    ccwBorderEdges=0;
    kk=1;
    for ii=1:nRefine % Bottom side
        jj=1;
        linInd=ii+(jj-1)*(nRefine);
        ccwBorderEdges(kk)=linInd;
        kk=kk+1;
    end
    
    for ii=1:nRefine % Right side
        jj=nRefine+1;
        linInd=ii+(jj-1)*(nRefine)+(nRefine*(nRefine+1));
        ccwBorderEdges(kk)=linInd;
        kk=kk+1;
    end
    for ii=nRefine:-1:1 % Top side
        jj=nRefine+1;
        linInd=ii+(jj-1)*(nRefine);
        ccwBorderEdges(kk)=linInd;
        kk=kk+1;
    end
    
    for ii=nRefine:-1:1 % Left side
        jj=1;
        linInd=ii+(jj-1)*(nRefine)+(nRefine*(nRefine+1));
        ccwBorderEdges(kk)=linInd;
        kk=kk+1;
    end
    
end

%% Base Structural Operation

function [edge]=AddEdgeStructure(index,cellIndex,vertexIndex,orientation)
    % Creates the strucure for an edge
    
    edge.index=index;
    edge.cellindex=cellIndex;
    edge.vertexindex=vertexIndex;
    edge.orientation=orientation;
%     edge.boundaryisHalf=false;
%     edge.boundaryis0=false;
%     edge.boundaryis1=false;
%     edge.solidisIn0=false;
%     edge.solidnotIn0=false;
%     edge.solidisIn1=false;
    
end

function [cell]=AddCellStructure(index,fill)
    % creates the structure for a cell
    cell.index=index;
    cell.fill=fill;
%     cell.coord=[];
%     cell.vertexindex=[];
    cell.refineLvl=1;
    cell.refineVec=1;
    cell.isactive=1;
end

function [vertex]=AddVertexStructure(index,coord)
    % creates the structure for a vertex
    vertex.index=index;
    vertex.coord=coord;
    
end

%% Refine Fully triangular mesh

function [refGrid,origconnec]=RefineTriangular(baseGrid,nRefine)
    
    tempGrid=baseGrid;
    delField={'boundaryisHalf','boundaryis0','boundaryis1','solidisIn0',...
        'solidnotIn0','solidisIn1','fill'};
    tempGrid.edge=rmfield(tempGrid.edge,delField);
    origconnec.edge=struct([]);
    [origconnec.cell(1:numel(baseGrid.cell)).old]=deal(baseGrid.cell.index);
    [origconnec.cell(1:numel(baseGrid.cell)).new]=deal(baseGrid.cell.index);
    
    [tempGrid,connecstructinfo]=RefineCellViaCentre(tempGrid);
    newOld=[connecstructinfo.cell.old];
    for jj=1:numel(origconnec.cell)
        origconnec.cell(jj).new=[connecstructinfo.cell(FindObjNum([],...
            origconnec.cell(jj).new,newOld)).new];
    end
    
    for ii=1:nRefine
        
        [tempGrid,cellVert]=BreakEdgesHalf(tempGrid);
        [tempGrid,connecstructinfo]=AddNewCellsAndEdges(tempGrid,cellVert);
        
        newOld=[connecstructinfo.cell.old];
        for jj=1:numel(origconnec.cell)
            origconnec.cell(jj).new=[connecstructinfo.cell(FindObjNum([],...
                origconnec.cell(jj).new,newOld)).new];
        end
        
    end
    
    refGrid=tempGrid;
end

function [refGrid,connec]=RefineCellViaCentre(baseGrid)
    
    cellCentreGrid=CellCentredGrid(baseGrid);
    [edgeInd]=[baseGrid.edge.index];
    cellInd=[baseGrid.cell.index];
    vertInd=[baseGrid.vertex.index];
    
    maxEdgeInd=max(edgeInd)+1;
    maxCellInd=max(cellInd)+1;
    maxVertInd=max(vertInd)+1;
    
    connec.edge=struct([]);
    connec.cell=repmat(struct('old',[],'new',[]),size(baseGrid.cell));
    
    for ii=1:numel(cellCentreGrid)
        nVert=numel(baseGrid.vertex);
        nCell=numel(baseGrid.cell);
        nEdge=numel(baseGrid.edge);
        
        newCoord=mean(vertcat(cellCentreGrid(ii).vertex.coord));
        baseGrid.vertex(nVert+1).index=maxVertInd;
        baseGrid.vertex(nVert+1).coord=newCoord;
        
        cellEdgeInd=[cellCentreGrid(ii).edge.index];
        [cellOrderedVertex,cellOrderedEdges]=OrderBlockEdges(...
            vertcat(cellCentreGrid(ii).edge.vertexindex),...
            cellEdgeInd);
        if numel(cellOrderedEdges)>1
            error('Unknown error in orderblockedges')
        end
        nNewEdge=numel(cellCentreGrid(ii).edge);
        newCellInd=[cellCentreGrid(ii).index,maxCellInd:maxCellInd+1];
        [baseGrid.cell(nCell+1:nCell+nNewEdge-1)]=deal(baseGrid.cell(ii));
        for jj=1:nNewEdge-1
            baseGrid.cell(nCell+jj).index=newCellInd(jj+1);
            maxCellInd=maxCellInd+1;
        end
        connec.cell(ii).old=baseGrid.cell(ii).index;
        connec.cell(ii).new=newCellInd;
        
        edgeSub=FindObjNum([],cellEdgeInd(cellOrderedEdges{1}),edgeInd);
        for jj=1:nNewEdge
            logCell=baseGrid.edge(edgeSub(jj)).cellindex==newCellInd(1);
            baseGrid.edge(edgeSub(jj)).cellindex(logCell)=newCellInd(jj);
            baseGrid.edge(nEdge+jj)=AddEdgeStructure(maxEdgeInd,...
                newCellInd([mod(jj-2,nNewEdge)+1,jj]),...
                [cellOrderedVertex{1}(jj,1),maxVertInd],0.5);
            maxEdgeInd=maxEdgeInd+1;
        end
        maxVertInd=maxVertInd+1;
    end
    
    refGrid=baseGrid;
end


function [baseGrid,cellVert]=BreakEdgesHalf(baseGrid)
    
    vertexInd=[baseGrid.vertex.index];
    vertexCoord=vertcat(baseGrid.vertex.coord);
    cellInd=[0,baseGrid.cell.index];
    
    maxVertInd=max(vertexInd)+1;
    maxEdgeInd=max([baseGrid.edge.index])+1;
    cellVert.cell=repmat(struct('index',[],'vertex',[]),[numel(baseGrid.cell)+1,1]);
    cellVert.vertex=repmat(struct('index',[],'edge',[]),[numel(baseGrid.edge),1]);
    
    for ii=1:numel(baseGrid.edge)
        
        vertSub=FindObjNum([],baseGrid.edge(ii).vertexindex,vertexInd);
        cellSub=FindObjNum([],baseGrid.edge(ii).cellindex,cellInd);
        
        kk=numel(baseGrid.edge)+1;
        baseGrid.edge(kk)=baseGrid.edge(ii);
        baseGrid.edge(kk).index=maxEdgeInd;
        
        kkv=numel(baseGrid.vertex)+1;
        baseGrid.vertex(kkv).index=maxVertInd;
        
        baseGrid.vertex(kkv).coord=sum(vertexCoord(vertSub,:))/2;
        
         baseGrid.edge(kk).vertexindex(1)=maxVertInd;
         baseGrid.edge(ii).vertexindex(2)=maxVertInd;
         for jj=1:numel(cellSub)
             cellVert.cell(cellSub(jj)).index=cellInd(cellSub(jj));
             cellVert.cell(cellSub(jj)).vertex=[cellVert.cell(cellSub(jj)).vertex,maxVertInd];
         end
         cellVert.vertex(ii).index=maxVertInd;
         cellVert.vertex(ii).edge=[baseGrid.edge(ii).index,maxEdgeInd];
         
         maxVertInd=maxVertInd+1;
         maxEdgeInd=maxEdgeInd+1;
    end

end

function [refGrid,connec]=AddNewCellsAndEdges(baseGrid,cellVert)
    % Refines a triangular cell by linkin all the new vertices
    cellInd=[0,baseGrid.cell.index];
    edgeInd=[baseGrid.edge.index];
    
    maxEdgeInd=max(edgeInd)+1;
    maxCellInd=max(cellInd)+1;
    
    vertCellVert=[cellVert.vertex.index];
    vertCellCell=[cellVert.cell.index];
    nEdge=numel(baseGrid.edge);
    
    baseGrid.edge(end+3*numel(baseGrid.cell)).index=[];
    
    connec.edge=struct([]);
    connec.cell=repmat(struct('old',[],'new',[]),size(baseGrid.cell));
    
    for ii=1:numel(baseGrid.cell)
        
        subCellVert=FindObjNum([],baseGrid.cell(ii).index,vertCellCell);
        newCellInd=maxCellInd:maxCellInd+2;
        
        newEdgeInd=[maxEdgeInd:maxEdgeInd+2]';
        edgeCellIndNew=[ones(3,1)*baseGrid.cell(ii).index,newCellInd'];
        edgeVertIndNew=cellVert.cell(subCellVert).vertex([1 2;2 3;3 1]);
        
        for jj=1:3;
            nEdge=nEdge+1;
            baseGrid.edge(nEdge).index=newEdgeInd(jj);
            baseGrid.edge(nEdge).vertexindex=edgeVertIndNew(jj,:);
            baseGrid.edge(nEdge).cellindex=edgeCellIndNew(jj,:);
            baseGrid.edge(nEdge).orientation=0.5;
        end
        
        % for each new Edge need to find the two old edge which closes the
        % containment and create the new cell
        % these are the ones which share a vertex
        
        %edgeCellSub=ceil(FindObjNum([],baseGrid.cell(ii).index,cellEdgeInd)/2);
        for jj=1:3
            vertSub=FindObjNum([],edgeVertIndNew(jj,:),vertCellVert);
            edgeSub=FindObjNum([],[cellVert.vertex(vertSub).edge],edgeInd);
            vertEdgeSub=[baseGrid.edge(edgeSub).vertexindex];
            vertEdgeSub=vertEdgeSub((vertEdgeSub~=edgeVertIndNew(jj,1)) & ...
                (vertEdgeSub~=edgeVertIndNew(jj,2)));
            cellSim=FindIdenticalVector(vertEdgeSub');
            edgeSub=edgeSub(cellSim{cellfun(@numel,cellSim)>1});
            for kk=1:2
                baseGrid.edge(edgeSub(kk)).cellindex(...
                    baseGrid.edge(edgeSub(kk)).cellindex==...
                    baseGrid.cell(ii).index)=newCellInd(jj);
            end
            baseGrid.cell(end+1)=baseGrid.cell(ii);
            baseGrid.cell(end).index=newCellInd(jj);
        end
        connec.cell(ii).old=baseGrid.cell(ii).index;
        connec.cell(ii).new=[baseGrid.cell(ii).index,newCellInd];
        % And then create the new cells
        
        maxCellInd=maxCellInd+3;
        maxEdgeInd=maxEdgeInd+3;
    end
    refGrid=baseGrid;
end


%% Plot Functions
function []=CheckResultsRefineGrid(gridreshape)
    domainBounds=[-1 1 ; -1 1];
    
    [unstructured]=ModifReshape(gridreshape);
    
    
    figh=figure;
    axh=axes;
    hold on
    
    colString='bgcmyk';
    
    isEdgeSub=find(unstructured.edge.boundaryis0);
    for ii=1:length(isEdgeSub)
        PlotEdge(figh,axh,unstructured,isEdgeSub(ii),'b-')
    end
    
    isEdgeSub=find(~unstructured.edge.boundaryis0);
    for ii=1:length(isEdgeSub)
        PlotEdge(figh,axh,unstructured,isEdgeSub(ii),'b-')
    end
    
    
    axis equal
    axis([domainBounds(1,1:2) domainBounds(2,1:2)])
    
    
end

function []=PlotEdge(figh,axh,unstructured,subEdge,format)
    figure(figh)
    %axes(axh)
    
    vertices=unstructured.edge.vertexindex(subEdge,:);
    vertsub(1)=find(unstructured.vertex.index==vertices(1));
    vertsub(2)=find(unstructured.vertex.index==vertices(2));
    coord=unstructured.vertex.coord(vertsub,:);
    
    plot(coord(:,1),coord(:,2),format)
    text(mean(coord(:,1)),mean(coord(:,2)),int2str(unstructured.edge.index(subEdge)),'color','b')
    text(mean(coord(1,1)),mean(coord(1,2)),int2str(vertices(1)),'color','g')
    text(mean(coord(2,1)),mean(coord(2,2)),int2str(vertices(2)),'color','g')
end

function []=PlotLoop(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).vertex.coord;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotSubDiv(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)
    
    
    coord=loop(indexLoop).subdivision;
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotCell(figh,axh,unstructured,indexCell,format)
    figure(figh)
    axes(axh)
    
    
    coord=unstructured.cell.coord(indexCell,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end


%% Test Code

function []=ExecuteTestCode()
    % wrapper function for test code 
    
    % Load test Data
    load([cd,'\TestData\TestDat_GridRefinement.mat'])
    
    
    testResult=TestFunctionsCall(standardInput);
    
    
    
    % Written Report
    runTimeError=find(testResult==1);
    executionError=find(testResult==-1);
    noError=find(testResult==0);
    fprintf('\n\n TEST CODE REPORT \n\n')
    fprintf('%i of %i tests executed without errors\n',numel(noError),numel(testResult))    
    fprintf('%i of %i tests did not produce the expected results\n',numel(runTimeError),numel(testResult))    
    fprintf('%i of %i tests could not be executed due to an error\n',numel(executionError),numel(testResult))
    fprintf('\n\n END OF REPORT \n\n')

end

function testResult=TestFunctionsCall(standardInput)
    % Function calling all the required Test functions
    testResult=[];
    kk=1;

    % Test Statements
    [testResult(kk)]=TestSplitEdge(standardInput);
    kk=kk+1;
    
    
end

function [testResult]=TestSplitEdge(standardInput)
    % Test Code Template
    % Expected Result
    expectres1(1)=AddEdgeStructure(25,[4,1],[5,17],[0,0]);
    expectres1(2)=AddEdgeStructure(26,[4,1],[17,6],[0,0]);
    expectres2=AddVertexStructure(17,[-2/3,-1/3]);
    % Function Call
    
    gridreshape=standardInput;
    edgeIndex=6;
    nRefine=2;
    
    try
        [newedges,newvertices]=SplitEdge(edgeIndex,nRefine,gridreshape);
    catch
        testResult=-1;
    end
    % Comparison
    testLog=isequal(expectres1,newedges)...
        && prod(abs([expectres2(:).coord]-[newvertices(:).coord])<1e-8)...
        && ([expectres2(:).index]==[newvertices(:).index]);
    if testLog
        testResult=0;
    else
        testResult=1;
    end
    
end

function [testResult,kk]=TestTemplate(kk,standardInput)
    % Test Code Template
    kk=kk+1;
    % Expected Result
    expectRes=0;
    % Function Call
    templateInput=standardInput;
    testResult=0;
    
    try
        outRes=Template(templateInput);
        % Comparison
    catch
        testResult=-1;
    end
    
    if testResult==0
        testLog=outRes==expectRes;
        if testLog
            testResult=0;
        else
            testResult=1;
        end
    end
end
