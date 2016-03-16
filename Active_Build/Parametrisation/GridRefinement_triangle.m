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



function [gridrefined,connectstructinfo,unstructuredrefined,loop]=...
        GridRefinement_triangle(gridreshape,nRefine,boundstr,typeRefine,execTest)
    % bridge function for the refinement function
    if ~exist('execTest','var'), execTest=true; end
    if ~exist('typeRefine','var'), typeRefine='grey'; end
    if execTest
        ExecuteTestCode();
    end
    
    [gridrefined,connectstructinfo]=RefineGrid(gridreshape,nRefine,typeRefine);
    
    [gridrefined]=EdgePropertiesReshape(gridrefined);
    isEdge=[gridrefined.edge(:).(boundstr{1})];
    cond=boundstr{3};
    [loop]=OrderSurfaceVertexReshape(gridrefined,isEdge,cond);
    
    [unstructuredrefined]=ModifReshape(gridrefined);
end

function [gridrefined,connecstructinfo]=RefineGrid(gridreshape,nRefine,typeRefine)
    % Function containing the grid refinement 
    
    [gridreshape.cell(:).coord]=deal([]);
    [gridreshape.cell(:).vertexindex]=deal([]);
    % Pre refinement processing
    [cellRefine]=IdentifyRefineCell(gridreshape,typeRefine);
    [edgeRefine]=IdentifyRefineEdges(gridreshape,cellRefine);
    [template,ccwtempinfo]=CreateCellRefinementTemplate(nRefine);
    
    % Refine Mesh
    [gridedgerefined,connecstructinfo]=...
        SplitSpecifiedEdges(gridreshape,edgeRefine,nRefine);
    
    
    
    [gridrefined,connecstructinfo.cell]=...
        SplitSpecifiedCells(gridedgerefined,cellRefine,template,ccwtempinfo);
%     CheckResultsRefineGrid(gridreshape)
%     CheckResultsRefineGrid(gridedgerefined)
%      CheckResultsRefineGrid(gridrefined)
    
end

%% Refinement preparation operations

function [cellRefine]=IdentifyRefineCell(gridreshape,typeRefine)
    % Identifies which cell need to be refined based on the type of
    % refinement call required
    % Returns the indices of the cells that must be refined
    
    switch typeRefine
        case 'all'
            cellRefine=[gridreshape.cell(:).index];
        case 'grey'
            cellRefineLog=[gridreshape.cell(:).fill]~=1 & [gridreshape.cell(:).fill]~=0;
            cellRefine=[gridreshape.cell(cellRefineLog).index];
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
    [templateTriangle,ccwinfoTriangle]=TriangularTemplate();
    cellIndices=[template.cell(:).index];
    [template,connecstructInfoTriangle]=...
        SplitSpecifiedCells(template,cellIndices,templateTriangle,ccwinfoTriangle);
    [edgeIndex]=[template.edge(:).index];
    [~,iOrder]=sort(edgeIndex);
    template.edge=template.edge(iOrder);
    [vertIndex]=[template.vertex(:).index];
    [~,iOrder]=sort(vertIndex);
    template.vertex=template.vertex(iOrder);
    cellIndex=[template.cell(:).index];
    newCellInd=cellIndex-(min(cellIndex)-1);
    cellConnec=vertcat(template.edge(:).cellindex);
    for ii=1:length(cellIndex)
        [xSub,ySub]=find(cellConnec==cellIndex(ii));
        template.cell(ii).index=newCellInd(ii);
        for jj=1:length(xSub)
            template.edge(xSub(jj)).cellindex(ySub(jj))=newCellInd(ii);
        end
        
    end
    
end%

function [gridreshape,connecstructinfo]=...
        SplitSpecifiedEdges(gridreshape,edgeIndices,nRefine)
    % Splits all the edges in the grid corresponding to a cell that will be
    % refined
    
    connecstructinfo.edge=[];
    
    nIndices=length(edgeIndices);
    
    for ii=1:nIndices
        
        [newedges,newvertices]=SplitEdge(edgeIndices(ii),nRefine,gridreshape);
        gridreshape.edge=[gridreshape.edge,newedges];
        gridreshape.vertex=[gridreshape.vertex,newvertices];
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
            vertexStream(ii:ii+1),baseedge.fill);
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
%     CheckResultsRefineGrid(adaptedtemplate)
    cellDelSub=FindObjNum(gridreshape.cell,cellInd);
    edgeDelSub=FindObjNum(gridreshape.edge,edgeList);
    edgeDelSub(edgeDelSub==0)=[];
    vertDelSub=FindObjNum(gridreshape.vertex,vertList);
    vertDelSub(vertDelSub==0)=[];
    
    
    gridreshape.cell(cellDelSub)=[];
    gridreshape.edge(edgeDelSub)=[];
    gridreshape.vertex(vertDelSub)=[];
    
    gridreshape.cell=[gridreshape.cell,adaptedtemplate.cell];
    gridreshape.vertex=[gridreshape.vertex,adaptedtemplate.vertex];
    gridreshape.edge=[gridreshape.edge,adaptedtemplate.edge];
    
    connectcellinfo.old=cellInd;
    connectcellinfo.new=cellList;
    
end

function [adaptedtemplate]=AdaptTemplate(template,edgeList,vertList,cellList,gridreshape,cellInd)
    % Adapts the template to match the specific indices of the case being
    % considered
    
    cornInd(1)=vertList(1);
    cornInd(2)=vertList(end);
    cornSub=FindObjNum(gridreshape.vertex,cornInd);
    baseCoord=vertcat(gridreshape.vertex(cornSub).coord);
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
    [cellIndTemplate]=[template.cell(:).index];
    
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
           
           edgetemplate(linInd)=AddEdgeStructure(indEdge,cellIndex,vertexIndex,[1,1]);
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
           
           edgetemplate(linInd)=AddEdgeStructure(indEdge,cellIndex,vertexIndex,[1,1]);
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

function [edge]=AddEdgeStructure(index,cellIndex,vertexIndex,fill)
    % Creates the strucure for an edge
    
    edge.index=index;
    edge.cellindex=cellIndex;
    edge.vertexindex=vertexIndex;
    edge.fill=fill;
    edge.boundaryisHalf=false;
    edge.boundaryis0=false;
    edge.boundaryis1=false;
    edge.solidisIn0=false;
    edge.solidnotIn0=false;
    edge.solidisIn1=false;
    
end

function [cell]=AddCellStructure(index,fill)
    % creates the structure for a cell
    cell.index=index;
    cell.fill=fill;
    cell.coord=[];
    cell.vertexindex=[];
end

function [vertex]=AddVertexStructure(index,coord)
    % creates the structure for a vertex
    vertex.index=index;
    vertex.coord=coord;
    
end

%% Triangular template
function [template,ccwinfo]=TriangularTemplate()
    [vertextemplate,ccwBorderVertices]=TemplateVerticesTriangle(1);
    [celltemplate]=TemplateCellTriangle();
    [edgetemplate,ccwBorderEdges]=TemplateEdgeTriangle();
    
    ccwinfo.edge=ccwBorderEdges;
    ccwinfo.cell=[];
    ccwinfo.vertex=ccwBorderVertices;
    
    template.vertex=vertextemplate;
    template.cell=celltemplate;
    template.edge=edgetemplate;
end

function [vertextemplate,ccwBorderVertices]=TemplateVerticesTriangle(nRefine)
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

function [celltemplate]=TemplateCellTriangle()
    % Creates the vertex template and returns the border
    
    for ii=1:2
       
            linInd=ii;
            [celltemplate(linInd)]=AddCellStructure(linInd,1);
            
        
    end 
end

function [edgetemplate,ccwBorderEdges]=TemplateEdgeTriangle()
    % Creates the edge template and returns the border
    
    cellIndex=[1,0;
        0,2;
        2,1;
        1,0;
        0,2];
    vertexIndex=[1,2;
        1,3;
        1,4;
        2,4;
        3,4];
    
    for ii=1:5
        
        edgetemplate(ii)=AddEdgeStructure(ii,cellIndex(ii,:),vertexIndex(ii,:),[1,1]);
    end
    [ccwBorderEdges]=[1,4,5,2];
    
end



%% Various Utilities

function []=template()
    
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
