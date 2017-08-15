function []=NURBSEngine(res)
    
    [nurbstruct,fullsnax]=ComputeRSVSNURBS(res.restartsnak.snaxel,res.grid);
    [loop]=OrderSurfaceSnaxel(fullsnax,fullsnax);
    u=linspace(0,1,3001);
    
    for ii=1:numel(nurbstruct)
        C=PlotNURBS(u,nurbstruct(ii).P,nurbstruct(ii).U,nurbstruct(ii).w,2);
        loop(ii).nurbs.pts=C;
        loop(ii).nurbs.ctrl=nurbstruct(ii).P;
    end
    
    figure,
        subplot(2,1,1)
        hold on
        subplot(2,1,2)
        hold on
    plotPoints= @(points,format) plot(points([1:end],1),points([1:end],2),format);
    for ii=1:numel(loop)
        subplot(2,1,1)
        plotPoints(loop(ii).nurbs.ctrl,'s-')
        plotPoints(loop(ii).nurbs.pts,'-')
        plotPoints(loop(ii).snaxel.coord,'*')
        subplot(2,1,2)
        [~,ptsNorm]=NormalDistance(loop(ii).snaxel.coord,loop(ii).nurbs.pts);
        plotPoints(ptsNorm,'*-')
    end
    
end

function [nurbstruct,fullsnax]=ComputeRSVSNURBS(snaxel,grid)
    % Translates a Snake into its NURBS representation
    % grid must contain fields, base, refined and connec
    nurbstruct=NURBSStructConstructor(2,0);
    
    [cellCentredCoarse]=CellCentredSnaxelInfoFromGrid(snaxel,grid);
    
    [snaxExtend,fullsnax]=BuildCoarseSnake(cellCentredCoarse);
    
    [gridcoarsen,coarsenconnec]=CoarsenGrid(grid);
    
    snakposition=PositionSnakes(snaxExtend,gridcoarsen);
    [snakposition]=SnaxelNormal2(snaxExtend,snakposition);
    [volumefraction,coeffstruct,cellCentredGrid]=VolumeFraction(snaxExtend,...
        snakposition,gridcoarsen,coarsenconnec,...
        CellCentreGridInformation(gridcoarsen),false(size(gridcoarsen.edge)));
    
    sub=FindObjNum([],[snaxExtend.snaxnext],[snaxExtend.index]);
    for ii=1:numel(snaxExtend)
        [nurbstructPart(ii)]=SnaxToNURBSPatch(snaxExtend([ii,sub(ii)]));
    end
    
   [nurbstruct]=BuildNurbsLoops(nurbstructPart,sub,2);
    
end

function [snaxExtend,fullsnax]=BuildCoarseSnake(cellCentredCoarse)
    % Extracts only the coarse Snaxel
    
    snaxExtend=[cellCentredCoarse(:).snaxel];
    [snaxExtendIndex,uniqSnax]=unique([snaxExtend.index]);
    snaxExtend=snaxExtend(uniqSnax);
    
    [snaxExtend]=CalculateSnaxelTangent(snaxExtend);
    [snaxExtend(:).cell]=deal([]);
    for ii=1:numel(cellCentredCoarse)
        if ~isempty(cellCentredCoarse(ii).snaxel)
            sub=FindObjNum([],[cellCentredCoarse(ii).snaxel.index],snaxExtendIndex);
            for jj=1:numel(sub)
                snaxExtend(sub(jj)).cell=[snaxExtend(sub(jj)).cell,cellCentredCoarse(ii).index];
            end
        end
    end
    fullsnax=snaxExtend;
    delIndex=snaxExtendIndex(~logical([snaxExtend.isborder]));
    snaxExtend=DeleteSnaxel(snaxExtend,delIndex);
end

function [nurbstruct]=BuildNurbsLoops(nurbstructPart,nextsub,deg)
    
     % loop the nurbs patches together
    nLoops=1;
    indList=1:numel(nurbstructPart);
    while any(indList)
        nurbstruct(nLoops)=NURBSStructConstructor(deg,0);
        currii=find(indList~=0,1,'first');
        flag=true;
        while flag
            nurbstruct(nLoops).P=[nurbstruct(nLoops).P;nurbstructPart(currii).P(1:2,:)];
            nurbstruct(nLoops).w=[nurbstruct(nLoops).w;nurbstructPart(currii).w(1:2)];
            indList(indList==currii)=0;
            currprec=currii;
            currii=nextsub(currii);
            flag=any(currii==indList);
        end
        nurbstruct(nLoops).P=[nurbstruct(nLoops).P;nurbstructPart(currprec).P(3,:)];
        nurbstruct(nLoops).w=[nurbstruct(nLoops).w;nurbstructPart(currprec).w(3)];
        d=nurbstruct(nLoops).d;
        nurbstruct(nLoops).U=[zeros(d,1);(ceil((0:size(nurbstruct(nLoops).P,1))/...
            2)/ceil(size(nurbstruct(nLoops).P,1))/2)';ones(d,1)];
        nLoops=nLoops+1;
    end
end

function nurbstruct=NURBSStructConstructor(deg,npts)
    
    nurbstruct.P=zeros(npts,2);
    nurbstruct.U=zeros((npts+1+deg)*sign(npts),1);
    nurbstruct.w=ones(npts,1);
    nurbstruct.d=deg;
    
end

function [vArc,P1,r,w,alpha,beta]=ComputeArcArea(P0,P2,a)
    % accepts row vectors
    normVec=@(v) sqrt(sum(v.^2));
    N=(P2-P0)*[0 -1;1 0];
    N0=P0+(P2-P0)/2;
    P1=N0+a*N;
    
    w=normVec(N0-P0)/normVec(P0-P1);
    alpha=asin(normVec(N0-P0)/normVec(P0-P1));
    beta=pi/2-alpha;
    r=normVec(N0-P0)/cos(alpha);
    
    vArc=(beta*normVec(P0-P1)^2/normVec(P1-N0)...
        -normVec(P0-N0))*normVec(P0-N0)^2/normVec(P1-N0);
    vArc(isnan(vArc))=0;
    vArc=sign(a)*vArc;
    
end

function [cellCentredCoarse]=CellCentredSnaxelInfoFromGrid(snaxel,grid)
    
    cellCentredBase=CellCentreGridInformation(grid.base);
    cellCentredRef=CellCentreGridInformation(grid.refined);
    
     [cellCentredCoarse]=CellCentredSnaxelInfo(snaxel,grid.refined,...
        cellCentredRef,cellCentredBase,grid.connec);
end

function [nurbstruct]=SnaxToNURBSPatch(snaxel)
    nurbstruct=NURBSStructConstructor(2,3);
    
    intersectFunc=@(c1,c2,v1,v2) ([v1',-v2']\(-c1'+c2'))';
    %ab=([snaxel(1).tangent',-snaxel(2).tangent']\(-snaxel(1).coord'+snaxel(2).coord'))';
    
    Pm=(snaxel(1).coord+snaxel(2).coord)/2;
    vecNorm=([0 1;-1 0]*(Pm-snaxel(1).coord)')';
    
    ab=intersectFunc(snaxel(1).coord,Pm,snaxel(1).tangent,vecNorm);
    Pest(1,1:2)=snaxel(1).tangent*ab(1)+snaxel(1).coord;
    ab=intersectFunc(snaxel(2).coord,Pm,snaxel(2).tangent,vecNorm);
    Pest(2,1:2)=snaxel(2).tangent*ab(1)+snaxel(2).coord;
    
    curv(1)=sqrt(sum(snaxel(1).curv.^2));
    curv(2)=sqrt(sum(snaxel(2).curv.^2));
    curv=repmat(reshape(curv,[2,1]),[1,2]);
    
    nurbstruct.P(2,:)=sum(Pest./curv)./(sum(1./curv));
%     nurbstruct.P(2,:)=(Pest1+Pest2)/2;
    nurbstruct.P(1,:)=snaxel(1).coord;
    nurbstruct.P(3,:)=snaxel(2).coord;
    
    
    d(1)=sqrt(sum((nurbstruct.P(2,:)-nurbstruct.P(1,:)).^2));
    d(2)=sqrt(sum((nurbstruct.P(2,:)-nurbstruct.P(3,:)).^2));
    d(3)=sqrt(sum((nurbstruct.P(3,:)-nurbstruct.P(1,:)).^2));
    
    w(1)=d(3)/2/d(1);
    w(2)=d(3)/2/d(2);
    
    nurbstruct.w(2)=mean(w);
    
end
%% Classic NURBS problems
function []=ComputeCircularNURBS()
    
    
    knots=[0 0 0 0.45 0.5 0.95 1 1 1];
    pts=[0 0; 0.2 0; 0.5 0.5; 0.8 0;1 0];
    w=[1 1 1 1 1]';
    deg=2;
    plotPoints= @(points,form) plot(points([1:end],1),points([1:end],2),form);
    u=linspace(0,1,501);
    
    hold on;
    %plotPoints(pts,'s-');
    [C]=PlotNURBS(u,pts,knots,w,deg);
    
    %plotPoints(C,'-')
    
    % circle from NURBS
    deg=2;
    knots=[0 0 0 1/4 1/4 4/8 1/2 3/4 3/4 1 1 1 ];
    pts=[ 1 0.5; 1 0; 0.5 0 ; 0 0; 0 0.5 ; 0 1; 0.5 1 ; 1 1; 1 0.5];
    a=1/sqrt(2);
    w=[1 a 1 a 1 a 1 a 1]';
    [C]=PlotNURBS(u,pts,knots,w,deg);
    figure,
    subplot(1,2,1)
    hold on
    plotPoints(pts,'s-');
    plotPoints(C,'-')
    axis equal
    
    % circle from NURBS
    deg=2;
    knots=[0 0 0 1/4 1/4 4/8 1/2 3/4 3/4 1 1 1 ];
    pts=[ 1 0.5; 1 0; 0.5 0 ; 0 0; 0 0.5 ; 0 1; 0.5 1 ; 1.5 1; 1.5 2];
    a=1/sqrt(2);
    w=[1 a 1 a 1 a 1 a 1]';
    [C]=PlotNURBS(u,pts,knots,w,deg);
    
    subplot(1,2,2), hold on
    plotPoints(pts,'s-');
    plotPoints(C,'-')
    axis equal
    
end

function [C]=PlotNURBS(u,pts,knots,w,deg)
    
    Z=zeros(size(pts,1),numel(u));
    for ii=1:size(pts,1)
        Z(ii,:) = SplineBasis(ii-1,deg,u,knots);
    end
    
    % One shot calculation of NURBS
    C=(Z'*diag(w)*pts)./repmat((Z'*w),[1,size(pts,2)]);
    
end

function Z = SplineBasis(i,p,u,U)  % from Dom
    % recursive definition of the polynomial bases used for B-Splines and
    % NURBS
    % (c) Dominic Masters - 2017
    
    i = i+1;
    if p==0
        % Z=zeros(1,length(u));
        for j=1:length(u);
            if j==length(u)
                if u(j)>=U(i) && u(j)<=U(i+1)
                    Z(j)=1;
                else
                    Z(j) = 0;
                end
            else
                if u(j)>=U(i) && u(j)<U(i+1)
                    Z(j)=1;
                else
                    Z(j)=0;
                end
            end
        end
    else
        %         Z1=((u-U(i))/(U(i+p)-U(i))).*N(i-1,p-1,u,U);
        %         Z2=((U(i+p+1)-u)/(U(i+p+1)-U(i+1))).*N((i+1)-1,p-1,u,U);
        %         if isnan(Z1)==1;
        %             if isnan(Z2)==1;
        %                 Z=zeros(size(u));
        %             else
        %                 Z=Z2;
        %             end
        %         else
        %             if isnan(Z2)==1;
        %                 Z=Z1;
        %             else
        %                 Z=Z1+Z2;
        %             end
        %         end
        
        t1n = (u-U(i));
        t1d = (U(i+p)-U(i));
        t2n = (U(i+p+1)-u);
        t2d = (U(i+p+1)-U(i+1));
        
        b1 = SplineBasis(i-1,p-1,u,U);
        b2 = SplineBasis((i+1)-1,p-1,u,U);
        
        if t1d ~= 0
            t1 = t1n/t1d;
        else
            t1 = 0;
        end
        
        if t2d ~= 0
            t2 = t2n/t2d;
        else
            t2 = 0;
        end
        
        Z = t1 .* b1 + t2 .* b2;
        
        
        
    end
end

%% from Snakes


function [singlesnaxel]=ModifyConnection(singlesnaxel,connecRemove,connecReplace)
    % Removes a connection to a snaxel from another snaxel and replaces
    % it with the replacement connection
    refresh=singlesnaxel.connectivity==connecRemove;
    isNext=singlesnaxel.snaxnext==connecRemove;
    isPrec=singlesnaxel.snaxprec==connecRemove;
    if sum(refresh)==0
        error('invalid connection to break');
    end
    if isNext && isPrec
        warning('Topology Collapsing');
    end
    % then replace that number to the snaxel that is being introduced
    singlesnaxel.connectivity(refresh)=connecReplace;
    if isNext
        singlesnaxel.snaxnext=connecReplace;
    elseif isPrec
        singlesnaxel.snaxprec=connecReplace;
    end
end

function snaxel=DeleteSnaxel(snaxel,delIndex)
    
    numDel=length(delIndex);
    snaxelIndices=[snaxel(:).index];
    delSnaxSub=zeros([1,numDel]);
    for ii=1:numDel
        delSnaxSub(ii)=FindObjNum(snaxel,delIndex(ii),snaxelIndices);
        
        connect=snaxel(delSnaxSub(ii)).connectivity;
        targSnaxSub=FindObjNum(snaxel,connect,snaxelIndices);
        
        
        for jj=1:length(targSnaxSub)
            singlesnaxel=snaxel(targSnaxSub(jj));
            connecReplace=connect(connect~=connect(jj));
            connecRemove=delIndex(ii);
            if ~isempty(connecReplace)
                [singlesnaxel]=ModifyConnection(singlesnaxel,connecRemove,connecReplace);
                snaxel(targSnaxSub(jj))=singlesnaxel;
            else
                disp('Topology Collapse')
            end
            
        end
        
    end
    
    snaxel(delSnaxSub)=[];
    
end


%% Contour Normal Calculation

function [snakposition]=SnaxelNormal2(snaxel,snakposition)
    % Calculates the normal at the Snaxel (According to snakes with topology control)
    
    
    [contourStruct]=ContourNormal2(snaxel,snakposition);
    nSnax=length(snaxel);
    
    ind1=[contourStruct(:).index1];
    ind2=[contourStruct(:).index2];
    
    for ii=1:nSnax
        contourVecSub1=FindObjNum([],snaxel(ii).index,ind1);
        contourVecSub2=FindObjNum([],snaxel(ii).index,ind2);
        contourVecSub=[contourVecSub1,contourVecSub2];
        contourVecSub(contourVecSub==0)=[];
        
        snakposition(ii).normvector=...
            vertcat(contourStruct(contourVecSub).vector);
        
        indexMat=[[contourStruct(contourVecSub).index1]',...
            [contourStruct(contourVecSub).index2]'];
        isPrecConn=logical(sum(snaxel(ii).snaxprec==indexMat,2));
        
        isNextConn=logical(sum(snaxel(ii).snaxnext==indexMat,2));
        
        snakposition(ii).vectornext=contourStruct(contourVecSub(isNextConn)).vector;
        snakposition(ii).vectorprec=contourStruct(contourVecSub(isPrecConn)).vector;
        
    end
    
end

function [contourStruct]=ContourNormal2(snaxel,snakposition)
    % Calculates the normal of the contour
    
    snaxIndPos=[snakposition(:).index];
    [contourStruct]=ExtractContourSnaxel(snaxel);
    % Normal Case
    for ii=1:length(contourStruct)
        indices=[contourStruct(ii).index1,contourStruct(ii).index2];
        
        [contourStruct(ii).vector]=ContourEdgeNormal(snakposition,...
            indices,snaxIndPos);
    end
    % Backup Method
    ind1=[contourStruct(:).index1];
    ind2=[contourStruct(:).index2];
    for ii=1:length(contourStruct)
        if sum(abs(contourStruct(ii).vector))==0
            
            [contourStruct(ii).vector]=NormalContourAlternateMethods(snakposition,...
                ii,contourStruct,ind1,ind2,snaxIndPos);
        end
    end
    for ii=1:length(contourStruct)
        [contourStruct(ii).vector]=contourStruct(ii).vector/...
            norm(contourStruct(ii).vector);
    end
    
end

function [contourStruct]=ExtractContourSnaxel(snaxel)
    % Extracts the snaxel connected for each  contour edge
%     snaxInd=[snaxel(:).index]';
%     snaxConnects=vertcat(snaxel(:).connectivity);
%     snaxContour1=[snaxInd,snaxConnects(:,1)];
%     snaxContour2=[snaxInd,snaxConnects(:,2)];
%     snaxContour=[snaxContour1;snaxContour2];
%     cellSimilar=FindIdenticalVector(snaxContour);
%     
%     arraySimilar=vertcat(cellSimilar{:});
%     snaxContour=snaxContour(arraySimilar(:,1),:);
%     for ii=length(snaxContour(:,1)):-1:1
%         contourStruct(ii).index1=snaxContour(ii,1);
%         contourStruct(ii).index2=snaxContour(ii,2);
%     end

    [contourStruct(1:numel(snaxel)).index1]=deal(snaxel(:).index);
    [contourStruct(1:numel(snaxel)).index2]=deal(snaxel(:).snaxnext);

end

function [normalVector]=ContourEdgeNormal(snakposition,indices,snaxIndPos)
    % Calculates the normal of the contour
    
    snaxSubPos=FindObjNum([],indices,snaxIndPos);
    
    [normalVector]=NormalContourBaseMethods(snakposition,snaxSubPos);
    baseVectors=vertcat(snakposition(snaxSubPos).vector);
    [normalVector]=TestNormalOutPointing(baseVectors,normalVector);
    
end

function [normalVectorOut]=TestNormalOutPointing(baseVectors,normalVector)
    % Tests that the normal points out in provided a vector pointing out is
    % specified
    for jj=1:length(baseVectors(:,1))
        testDirNV(jj)=dot(baseVectors(jj,:),normalVector);
    end
    
    if sum(testDirNV<0)
        normalVectorOut=-normalVector;
    elseif sum(testDirNV==0)==2
        normalVectorOut=[0 0];
    else
        normalVectorOut=normalVector;
    end
    
end

function [normalVector]=NormalContourBaseMethods(snakposition,snaxSubPos)
    
    coord1=snakposition(snaxSubPos(1)).coord;
    coord2=snakposition(snaxSubPos(2)).coord;
    tanVec=coord1-coord2;
    % Normal Cases for calculation of
    if sum(abs(tanVec))~=0
        % if the tangential vector is not the 0 vector
        normalVector=CalcNormVec2DClockWise(tanVec);
    else
        % else use the sum of direction vectors
        normalVector=snakposition(snaxSubPos(1)).vector+...
            snakposition(snaxSubPos(2)).vector;
    end
    
    
end

function [normalVector]=NormalContourAlternateMethods(snakposition,currentCont,...
        contourStruct,ind1,ind2,snaxIndPos)
    % Alternate method relying on adjacent edge contours to define the
    % normal vector
    
    indices=[contourStruct(currentCont).index1,contourStruct(currentCont).index2];
    otherCont1=FindObjNum([],indices,ind1);
    otherCont2=FindObjNum([],indices,ind2);
    otherCont1(otherCont1==currentCont)=[];
    otherCont2(otherCont2==currentCont)=[];
    otherCont=[otherCont1',otherCont2'];
    otherCont(otherCont==0)=[];
    adjacentVectors=vertcat(contourStruct(otherCont).vector);
    outPointVec=sum(adjacentVectors);
    snaxSubPos=FindObjNum([],indices,snaxIndPos);
    
    coord1=snakposition(snaxSubPos(1)).coord;
    coord2=snakposition(snaxSubPos(2)).coord;
    tanVec=coord1-coord2;
    % Normal Cases for calculation of
    if sum(abs(tanVec))~=0
        % if the tangential vector is not the 0 vector
        normalVector=CalcNormVec2DClockWise(tanVec);
    else
        % else use the sum of direction vectors
        normalVector=CalcNormVec2DClockWise(snakposition(snaxSubPos(1)).vector);
    end
    [normalVector]=TestNormalOutPointing(outPointVec,normalVector);
    
end