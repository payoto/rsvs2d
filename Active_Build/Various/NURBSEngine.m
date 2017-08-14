function []=NURBSEngine(res0,resL)
    
    [nurbstruct]=ComputeRSVSNURBS(resL.restartsnak.snaxel,res0.grid);
    
    u=linspace(0,1,3001);
    loop=resL.loop;
    for ii=1:numel(nurbstruct)
        C=PlotNURBS(u,nurbstruct(ii).P,nurbstruct(ii).U,nurbstruct(ii).w,2);
        loop(ii).nurbs.pts=C;
        loop(ii).nurbs.ctrl=nurbstruct(ii).P;
    end
    
    figure,
    hold on
    plotPoints= @(points,format) plot(points([1:end],1),points([1:end],2),format);
    for ii=1:numel(loop)
        
        plotPoints(loop(ii).nurbs.ctrl,'bs-')
        plotPoints(loop(ii).nurbs.pts,'r-')
        plotPoints(loop(ii).snaxel.coord,'g*')
        
    end
    
end

function [nurbstruct]=ComputeRSVSNURBS(snaxel,grid)
    % Translates a Snake into its NURBS representation
    % grid must contain fields, base, refined and connec
    nurbstruct=NURBSStructConstructor(2,0);
    
    [cellCentredCoarse]=CellCentredSnaxelInfoFromGrid(snaxel,grid);
    
    
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
    delIndex=snaxExtendIndex(~logical([snaxExtend.isborder]));
    snaxExtend=DeleteSnaxel(snaxExtend,delIndex);
    
    sub=FindObjNum([],[snaxExtend.snaxnext],[snaxExtend.index]);
    for ii=1:numel(snaxExtend)
        [nurbstructPart(ii)]=SnaxToNURBSPatch(snaxExtend([ii,sub(ii)]));
    end
    
    % loop the nurbs patches together
    nLoops=1;
    indList=1:numel(nurbstructPart);
    while any(indList)
        nurbstruct(nLoops)=NURBSStructConstructor(2,0);
        currii=find(indList~=0,1,'first');
        flag=true;
        while flag
            nurbstruct(nLoops).P=[nurbstruct(nLoops).P;nurbstructPart(currii).P(1:2,:)];
            nurbstruct(nLoops).w=[nurbstruct(nLoops).w;nurbstructPart(currii).w(1:2)];
            indList(indList==currii)=0;
            currprec=currii;
            currii=sub(currii);
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
    Pest1=snaxel(1).tangent*ab(1)+snaxel(1).coord;
    ab=intersectFunc(snaxel(2).coord,Pm,snaxel(2).tangent,vecNorm);
    Pest2=snaxel(2).tangent*ab(1)+snaxel(2).coord;
    
    nurbstruct.P(2,:)=(Pest1+Pest2)/2;
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
        singlesnaxel.tovertex=connecReplace;
    elseif isPrec
        singlesnaxel.snaxprec=connecReplace;
        singlesnaxel.fromvertex=connecReplace;
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