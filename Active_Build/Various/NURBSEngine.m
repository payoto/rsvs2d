function [nurbstruct]=NURBSEngine(runType,snaxel,snaxgrid,varargin)
    % This function Generates and plots NURBS and corresponding snakes;
    
    switch runType
        case 'NURBSgen'
            [nurbstruct,fullsnax]=ComputeRSVSNURBS(snaxel,snaxgrid,'area');
            [nurbstruct(2,:),fullsnax]=ComputeRSVSNURBS(snaxel,snaxgrid,'grad');
            [loop]=OrderSurfaceSnaxel(fullsnax,fullsnax);
            loop=repmat(reshape(loop,[1,numel(loop)]),[2,1]);
            figName='Area and Gradient NURBS';
            cellStr={'area','grad'};
        case 'snakedense'
            [nurbstruct,fullsnax]=ComputeRSVSNURBS(snaxel{1},snaxgrid{1},'area');
            [loop]=OrderSurfaceSnaxel(fullsnax,fullsnax);
            loop=repmat(reshape(loop,[1,numel(loop)]),[numel(snaxel),1]);
            for ii=2:numel(snaxel)
                [nurbstruct(ii,:),fullsnax]=ComputeRSVSNURBS(snaxel{ii},snaxgrid{ii},'area');
                [loop(ii,:)]=OrderSurfaceSnaxel(fullsnax,fullsnax);
            end
            figName='Snake Density effect on NURBS';
            cellStr=varargin{1};
    end
    
    u=linspace(0,1,3001);
    for jj=1:size(nurbstruct,1)
        for ii=1:size(nurbstruct,2)
            C=PlotNURBS(u,nurbstruct(jj,ii).P,nurbstruct(jj,ii).U,nurbstruct(jj,ii).w,2);
            loop(jj,ii).nurbs.pts=C;
            loop(jj,ii).nurbs.ctrl=nurbstruct(jj,ii).P;
        end
    end
    
    h=figure('Name',figName,'Position',[200 200 1000 600]);
    ax(1)=subplot(2,1,1);
    hold on
    ax(2)=subplot(2,1,2);
    hold on
    plotPoints= @(points,format) plot(points([1:end],1),points([1:end],2),format);
    
    c=get(gca,'colororder');
    for ii=1:size(loop,2)
        subplot(2,1,1)
        kk=1;
        ckk=0;
        for jj=1:size(loop,1)
            subplot(2,1,1)
            ckk=1+ckk;
            l(kk)=plotPoints(loop(jj,ii).snaxel.coord,'+');
            l(kk).Color=c(mod(ckk-1,size(c,1))+1,:);
            l(kk).DisplayName='snaxel';kk=1+kk;
            ckk=1+ckk;
            l(kk)=plotPoints(loop(jj,ii).nurbs.ctrl,'s--');
            l(kk).Color=c(mod(ckk-1,size(c,1))+1,:);
            l(kk).DisplayName=['NURBS control - ',cellStr{jj}];kk=1+kk;
            l(kk)=plotPoints(loop(jj,ii).nurbs.pts,'-');
            l(kk).Color=c(mod(ckk-1,size(c,1))+1,:);
            l(kk).DisplayName=['NURBS - ',cellStr{jj}];kk=1+kk;
            subplot(2,1,2)
            [~,ptsNorm]=NormalDistance(loop(jj,ii).snaxel.coord,loop(jj,ii).nurbs.pts);
            l(kk)=plotPoints(ptsNorm,'*-');
            l(kk).Color=c(mod(ckk-1,size(c,1))+1,:);
            l(kk).DisplayName=['Distance - ',cellStr{jj}];kk=1+kk;
        end
    end
    ax(3)=legend(l);
    %ax=columnlegend(3,{l(:).DisplayName},'legend',legend_h,'object',object_h);
    
    ax(3).Orientation='horizontal';
    ax(3).Location='SouthOutside';
    ax(2).YScale='log';
    ax(1).OuterPosition=[0 0.35 1 0.65];
    ax(1).Position([1 3])=[ 0.05 0.90];
    ax(2).OuterPosition=[0 0.05 1 0.28];
    ax(2).Position([1 3])=[ 0.05 0.90];
    ax(2).YLim=[1e-8 1e-2];
    ax(2).YTick=10.^-[8:-2:1];
    ax(2).YGrid='on';
end

function []=PlotNurbStructLoop(loop)
    
end

function [nurbstruct,fullsnax]=ComputeRSVSNURBS(snaxel,grid,NURBScalc)
    % Translates a Snake into its NURBS representation
    % grid must contain fields, base, refined and connec
    nurbstruct=NURBSStructConstructor(2,0);
    
    if ~isfield(grid.connec.cell,'newCellInd')
        [grid.connec.cell.newCellInd]=deal(grid.connec.cell.new);
        [grid.connec.cell.oldCellInd]=deal(grid.connec.cell.old);
    end
    [cellCentredCoarse]=CellCentredSnaxelInfoFromGrid(snaxel,grid);
    
    [snaxExtend,fullsnax]=BuildCoarseSnake(cellCentredCoarse);
    sub=FindObjNum([],[snaxExtend.snaxnext],[snaxExtend.index]);
    
    
    switch NURBScalc
        case 'area'
            [gridcoarsen,coarsenconnec]=CoarsenGrid(grid);
            snakposition=PositionSnakes(snaxExtend,gridcoarsen);
            [snakposition]=SnaxelNormal2(snaxExtend,snakposition);
            [volumefraction,coeffstruct,cellCentredGrid]=VolumeFraction(snaxExtend,...
                snakposition,gridcoarsen,coarsenconnec,...
                CellCentreGridInformation(gridcoarsen),false(size(gridcoarsen.edge)));
            deltaVolume=[cellCentredGrid.fill].*[volumefraction.totalvolume]-...
                [volumefraction.totalfraction];
            nurbstructPart=GenerateNURBSArea(snaxExtend,deltaVolume,gridcoarsen);
        case 'grad'
            for ii=1:numel(snaxExtend)
                [nurbstructPart(ii)]=SnaxGradToNURBSPatch(snaxExtend([ii,sub(ii)]));
            end
            
        case 'errmin'
            % Build here a NURBS for error minimisation with the Snake
    end
    
    [nurbstruct]=BuildNurbsLoops(nurbstructPart,sub,2);
    
end

function [nurbpatch]=GenerateNURBSArea(snaxel,deltaVolume,refinedgrid)
    
    normMat=@(v) sqrt(sum((v(1,:)-v(2,:)).^2));
    nurbpatch=repmat(struct('snaxel',[0 0],'snaxsub',[0 0],'coord',[0 0 0 0],'cell',0,...
        'length',[0],'dvol',0,'P',zeros(3,2),'w',[0 0 0]'),size(snaxel));
    
    edgeInd=[refinedgrid.edge.index];
    sub=FindObjNum([],[snaxel.snaxnext],[snaxel.index]);
    
    deltaVolume=[deltaVolume];
    lengthScale=zeros(size(deltaVolume));
    
    for ii=1:numel(snaxel)
        nurbpatch(ii).snaxel=[snaxel(ii).index,snaxel(ii).snaxnext];
        nurbpatch(ii).snaxsub=[ii, sub(ii)];
        nurbpatch(ii).length=(normMat(vertcat(snaxel([ii, sub(ii)]).coord)));
        nurbpatch(ii).coord=[snaxel([ii, sub(ii)]).coord];
        
        intermCell=vertcat(refinedgrid.edge(FindObjNum([],...
            [snaxel([ii,sub(ii)]).edge],edgeInd)).cellindex);
        
        cellLog=(intermCell(1,:)'*ones(1,2)==ones(2,1)*intermCell(2,:));
        nurbpatch(ii).cell=intermCell(1,mod(find(cellLog)-1,2)+1);
    end
    kk=find(cellfun(@numel,{nurbpatch(:).cell})==1,1,'first');
    for ii=1:numel(nurbpatch)
        ind=mod(kk+ii-1,numel(nurbpatch))+1;
        prevind=mod(kk+ii-2,numel(nurbpatch))+1;
        if numel(nurbpatch(ind).cell)>1
            nurbpatch(ind).cell=nurbpatch(ind).cell(nurbpatch(ind).cell~=nurbpatch(prevind).cell);
        end
        lengthScale(nurbpatch(ind).cell)=lengthScale(nurbpatch(ind).cell)+nurbpatch(ind).length;
    end
    
    
    for ii=1:numel(snaxel)
        nurbpatch(ii).dvol=nurbpatch(ii).length*deltaVolume(nurbpatch(ii).cell)/lengthScale(nurbpatch(ii).cell);
        nurbpatch(ii).dvol(isnan(nurbpatch(ii).dvol))=0;
    end
    
    for ii=1:numel(nurbpatch)
        [nurbpatch(ii)]=MatchNURBSpatchArea(nurbpatch(ii));
    end
    
end

function [nurbpatch]=MatchNURBSpatchArea(nurbpatch)
    
    compAreaGS=@(a,P) -abs(ComputeArcArea(P(1:2),P(3:4),a)-P(5));
    
    b=[-1 1];
    bTest=[0 0];
    for ii=1:2
        bTest(ii)=ComputeArcArea(nurbpatch.coord(1:2),nurbpatch.coord(3:4),b(ii));
        while sign(b(ii))*bTest(ii)<sign(b(ii))*nurbpatch.dvol
            b(ii)=b(ii)*2;
            bTest(ii)=ComputeArcArea(nurbpatch.coord(1:2),nurbpatch.coord(3:4),b(ii));
        end
    end
    
    [np]=GoldenSection_func(b(1),b(2),1e-12,[nurbpatch.coord,nurbpatch.dvol],...
        compAreaGS);
    
    [vArc,P1,r,w,alpha,beta]=ComputeArcArea(nurbpatch.coord(1:2),nurbpatch.coord(3:4),np);
    
    nurbpatch.P=[nurbpatch.coord(1:2);P1;nurbpatch.coord(3:4)];
    nurbpatch.w=[1;w;1];
    
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

function [nurbstruct]=SnaxGradToNURBSPatch(snaxel)
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