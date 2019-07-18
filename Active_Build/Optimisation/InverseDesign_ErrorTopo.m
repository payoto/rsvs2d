%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%       Check inverse design match
%        for parametric snakes
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [] = InverseDesign_ErrorTopo()
%     FUNCTIONLIST allows local functions to be used globally once it has
%     been used.
%     
%     funcHandles=localfunctions;
%     funcDir=[cd,'\Automated_Function_Directory_ExecOptim'];
%     HeaderActivation(funcHandles,funcDir)
%     
% end

function [errorMeasure,h,targCoord,analysisLoop]=InverseDesign_ErrorTopo(paramoptim,loop, arg3)
    if nargin==2
        [errorMeasure,h,targCoord,analysisLoop]=ErrorTopo_ParamInput(paramoptim,loop);
    else
        analysisLoop=paramoptim;
        targCoord=loop;
        profileComp=arg3;
        [errorMeasure,h,targCoord,analysisLoop]=ErrorTopo(...
            analysisLoop,targCoord,profileComp);
    end
end

function [errorMeasure,h,targCoord,analysisLoop]=ErrorTopo_ParamInput(paramoptim,loop)
    
    varExtract={'aeroClass','aeroName','profileComp'};
    [aeroClass,aeroName,profileComp]=ExtractVariables(varExtract,paramoptim);
    varExtract={'typeLoop'};
    [typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    %profileComp='area';
    
    [analysisLoop,targPrep]=PrepareLoopCoord(loop,'area',typeLoop);
    
    switch aeroClass
        case 'NACA'
            [targCoord]=GenerateNacaCoordTOPO(targPrep(:,1),targPrep(:,2),aeroName);
        case 'lib'
            
            [targCoord]=GenerateLibCoord(targPrep(:,1),targPrep(:,2),aeroName);
        case 'loop'
            [targCoord]=LoadLoopData(aeroName);
        otherwise
            error('not coded yet')
    end
    [errorMeasure,h,targCoord,analysisLoop]=ErrorTopo(...
        analysisLoop,targCoord,profileComp);
end

function [errorMeasure,h,targCoord,analysisLoop]=ErrorTopo(...
        analysisLoop,targCoord,profileComp)
    
    switch profileComp
        case 'distance'
            error('Distance error measurement is not supported with topology')
        case 'area'
            [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaTopo(analysisLoop,targCoord);
        case 'areasquared'
            [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaSquaredTopo(analysisLoop,targCoord);
        case 'areadist'
            [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaRawDistTopo(analysisLoop,targCoord);
        case 'areapdist'
            [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaPRawDistTopo(analysisLoop,targCoord);
        case 'area2parea'
            [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaSquaredpAreaTopo(analysisLoop,targCoord);
        case 'normdist'
            [errorMeasure,modifiedDistance,analysisLoop]=CompareProfilesNormDistTopo(analysisLoop,targCoord);
            %ax.YScale='log';
        otherwise
            error('not coded yet')
    end
    plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
    h=figure;
    subplot(2,1,1)
    PlotLoop(analysisLoop,'coord')
    hold on
    PlotLoop(targCoord,'coord')
    if exist('indepLoop','var')
        patches=PlotLoop(indepLoop,'coord',1);
        [patches.FaceAlpha]=deal(0.2);
        [patches.LineStyle]=deal('none');
    end
    legend('snake points','Target points')
    ax=subplot(2,1,2);
    plotPoints(modifiedDistance)
    
end

function [errorMeasure,modifiedDistance,testLoop]=CompareProfilesNormDistTopo(testLoop,targLoop)
    % test intersection and internal with inpolygon
    % targLoop ->
    % testLoop |
    typeLoop='coord';
    
    [intersectTable]=BuildIntersectionTable(testLoop,targLoop,typeLoop);
    
    % once tested establish the appropriate objective function
    objChoice=all(any(intersectTable~=0,1));
    
    %
    if objChoice % use iterative area condition
        
        errorMeasure=cell([1,numel(testLoop)]);
        modifiedDistanceFull=cell([1,numel(testLoop)]);
        for ii=1:numel(testLoop)
            targActTest=find(intersectTable(ii,:));
            errorMeasureTemp=cell([1,numel(targActTest)]);
            for jj=1:numel(targActTest)
                analysisCoord=RemoveIdenticalConsecutivePoints(targLoop(...
                    targActTest(jj)).coord);
                testCoord=RemoveIdenticalConsecutivePoints(testLoop(ii).coord);
                [errorMeasureTemp{jj},modifiedDistanceTemp]=...
                    CompareProfilesDistance2(testCoord,analysisCoord);
                modifiedDistanceFull{ii}(:,1)=modifiedDistanceTemp(:,1);
                modifiedDistanceFull{ii}(:,1+jj)=modifiedDistanceTemp(:,2);
            end
            modifiedDistanceFull{ii}(:,2)=min(modifiedDistanceFull{ii}(:,2:end),[],2);
            modifiedDistanceFull{ii}(:,3:end)=[];
            errorMeasure{ii}.sum=sum(modifiedDistanceFull{ii}(:,2));
            errorMeasure{ii}.mean=mean(modifiedDistanceFull{ii}(:,2));
            errorMeasure{ii}.std=std(modifiedDistanceFull{ii}(:,2));
            errorMeasure{ii}.max=max(modifiedDistanceFull{ii}(:,2));
            errorMeasure{ii}.min=min(modifiedDistanceFull{ii}(:,2));
            testLoop(ii).localerror=modifiedDistanceFull{ii}(:,2);
        end
        for ii=2:numel(errorMeasure)
            errorMeasure{1}=opstruct('+',errorMeasure{[1,ii]});
        end
        errorMeasure=errorMeasure{1};
        modifiedDistance2=zeros([0 2]);
        for ii=1:numel(modifiedDistanceFull)
            modifiedDistance2=[modifiedDistance2;[1,0];modifiedDistanceFull{ii}];
        end
        modifiedDistance=modifiedDistance2;
%         [~,iSort]=sort(modifiedDistance(:,1));
%         modifiedDistance=modifiedDistance(iSort,:);
    else % use nearest neighbour aproach with area matching
        error('Profiles are not valid for distance matching')
    end
    
    
end

function [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaTopo(testLoop,targLoop)
    % test intersection and internal with inpolygon
    % targLoop ->
    % testLoop |
    typeLoop='coord';
    
    [intersectTable]=BuildIntersectionTable(testLoop,targLoop,typeLoop);
    
    % once tested establish the appropriate objective function
    objChoice=all(any(intersectTable~=0,1));
    
    %
    if objChoice % use iterative area condition
        [indepLoop]=AreaErrorTopo(testLoop,targLoop);
        
        targArea=0;
        for ii=1:numel(targLoop)
            targArea=targArea+abs(CalculatePolyArea(targLoop(ii).coord));
        end
        
        [errorMeasure,modifiedDistance,indepLoop]=IndepProfileError(indepLoop,targArea);
    else % use nearest neighbour aproach with area matching
        indepLoop=repmat(struct('coord',zeros([0 2])),[0 1]);
     [errA]=NotIntersectAreaCondition(targLoop,testLoop,typeLoop);
        [errorMeasure,modifiedDistance]=NotIntersectCondition(targLoop,testLoop,typeLoop,errA);
    end
    
    
end

function [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaSquaredTopo(testLoop,targLoop)
    % test intersection and internal with inpolygon
    % targLoop ->
    % testLoop |
    typeLoop='coord';
    
    [intersectTable]=BuildIntersectionTable(testLoop,targLoop,typeLoop);
    
    % once tested establish the appropriate objective function
    objChoice=all(any(intersectTable~=0,1));
    
    %
    if objChoice % use iterative area condition
        [indepLoop]=AreaErrorTopo(testLoop,targLoop);
        
        targArea=0;
        for ii=1:numel(targLoop)
            targArea=targArea+abs(CalculatePolyArea(targLoop(ii).coord));
        end
        funArea{1}=@(x) x.^2;
        funArea{2}=@(x) sqrt(x);
        [errorMeasure,modifiedDistance,indepLoop]=IndepProfileError(indepLoop,targArea,funArea);
    else % use nearest neighbour aproach with area matching
        indepLoop=repmat(struct('coord',zeros([0 2])),[0 1]);
        [errA]=NotIntersectAreaCondition(targLoop,testLoop,typeLoop);
        [errorMeasure,modifiedDistance]=NotIntersectCondition(targLoop,testLoop,typeLoop,errA);
    end
    
    
end

function [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaSquaredpAreaTopo(testLoop,targLoop)
    % test intersection and internal with inpolygon
    % targLoop ->
    % testLoop |
    typeLoop='coord';
    
    [intersectTable]=BuildIntersectionTable(testLoop,targLoop,typeLoop);
    
    % once tested establish the appropriate objective function
    objChoice=all(any(intersectTable~=0,1));
    
    %
    if objChoice % use iterative area condition
        [indepLoop]=AreaErrorTopo(testLoop,targLoop);
        
        targArea=0;
        for ii=1:numel(targLoop)
            targArea=targArea+abs(CalculatePolyArea(targLoop(ii).coord));
        end
        funArea{1}=@(x) x.^2;
        funArea{2}=@(x) sqrt(x);
        [errorMeasure,modifiedDistance,indepLoop]=IndepProfileError(indepLoop,targArea,funArea);
        [errorMeasure2,modifiedDistance2,indepLoop]=IndepProfileError(indepLoop,targArea);
        
        errorMeasure=opstruct(@(a,b)(a+b)/2,errorMeasure,errorMeasure2);
        modifiedDistance=(modifiedDistance+modifiedDistance2)/2;
    else % use nearest neighbour aproach with area matching
        indepLoop=repmat(struct('coord',zeros([0 2])),[0 1]);
        [errA]=NotIntersectAreaCondition(targLoop,testLoop,typeLoop);
        [errorMeasure,modifiedDistance]=NotIntersectCondition(targLoop,testLoop,typeLoop,errA);
    end
    
    
end

function [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaRawDistTopo(testLoop,targLoop)
    % test intersection and internal with inpolygon
    % targLoop ->
    % testLoop |
    typeLoop='coord';
    
    [intersectTable]=BuildIntersectionTable(testLoop,targLoop,typeLoop);
    
    % once tested establish the appropriate objective function
    objChoice=all(any(intersectTable~=0,1));
    
    targArea=0;
    for ii=1:numel(targLoop)
        targArea=targArea+abs(CalculatePolyArea(targLoop(ii).coord));
    end
    %
    if objChoice % use iterative area condition
        [indepLoop]=AreaErrorTopo(testLoop,targLoop);
        
        
        
        [errorMeasure,modifiedDistance,indepLoop]=IndepProfileError(indepLoop,targArea);
    else % use nearest neighbour aproach with area matching
        indepLoop=repmat(struct('coord',zeros([0 2])),[0 1]);
        
        [errorMeasure,modifiedDistance]=NotIntersectCondition(targLoop,testLoop,typeLoop,targArea);
    end
    
    
end

function [errorMeasure,modifiedDistance,indepLoop]=CompareProfilesAreaPRawDistTopo(testLoop,targLoop)
    % test intersection and internal with inpolygon
    % targLoop ->
    % testLoop |
    typeLoop='coord';
    
    [intersectTable]=BuildIntersectionTable(testLoop,targLoop,typeLoop);
    
    % once tested establish the appropriate objective function
    objChoice=all(any(intersectTable~=0,1));
    
    targArea=0;
    for ii=1:numel(targLoop)
        targArea=targArea+abs(CalculatePolyArea(targLoop(ii).coord));
    end
    %
    if objChoice % use iterative area condition
        [indepLoop]=AreaErrorTopo(testLoop,targLoop);
        
        
        
        [errorMeasure,modifiedDistance,indepLoop]=IndepProfileError(indepLoop,targArea);
        [errorMeasure2,modifiedDistance2]=NotIntersectCondition(targLoop,testLoop,typeLoop,0);
        errorMeasure=opstruct('+',errorMeasure,errorMeasure2);
    else % use nearest neighbour aproach with area matching
        indepLoop=repmat(struct('coord',zeros([0 2])),[0 1]);
        
        [errorMeasure,modifiedDistance]=NotIntersectCondition(targLoop,testLoop,typeLoop,targArea);
    end
    
    
end

function [out]=TestInternalOrIntersect(coord1,coord2)
    % takes in two coordinate lists:
    % returns 0 if they are not intersecting or iternal to each other
    % 1 if 1 is in 2
    % 2 if 2 is in 1
    % -1 if they intersect
    out=0;
    
    [in1,on1] = inpolygon(coord1(:,1),coord1(:,2),coord2(:,1),coord2(:,2));
    isIntersect= any(on1) || xor(any(in1),all(in1));
    if isIntersect
        out=-1;
    elseif any(in1)
        out=1;
    else
        [in2,~] = inpolygon(coord2(1,1),coord2(1,2),coord1(:,1),coord1(:,2));
        out=2*any(in2);
    end
end

function [errorMeasure,modifiedDistance]=NotIntersectCondition(targLoop,testLoop,typeLoop,errA)
    % NearestNeighbourCondition
    targCoord=vertcat(targLoop(:).coord);
    testCoord=vertcat(testLoop(:).(typeLoop));
    
    [idTest, dTarg] = knnsearch(testCoord,targCoord);
    
    modifiedDistance=[dTarg.^2];
    
    
    errorMeasure.sum=sum(modifiedDistance)+errA;
    errorMeasure.mean=mean(modifiedDistance)+errA;
    errorMeasure.std=std(modifiedDistance)+errA;
    errorMeasure.max=max(modifiedDistance)+errA;
    errorMeasure.min=min(modifiedDistance)+errA;
    modifiedDistance=[targCoord(:,1),modifiedDistance];
end

function [errA]=NotIntersectAreaCondition(targLoop,testLoop,typeLoop)
    % This term drives the profiles to an area similar as the target and
    % ensures an intersecting case will always be better in terms of
    % objective
    % err=2*Atarg+abs(Atarg-Atest)
    % This term tends to get the areas to be the same size and is always
    % larger than the maximum area error
    
    targA=zeros(size(targLoop));
    testA=zeros(size(testLoop));
    for ii=1:numel(testLoop)
        testA(ii)=abs(CalculatePolyArea(testLoop(ii).(typeLoop)));
    end
    for ii=1:numel(targLoop)
        targA(ii)=abs(CalculatePolyArea(targLoop(ii).coord));
    end
    
    errA=(2*sum(targA)+abs(sum(targA)-sum(testA)))/sum(targA);
end

function [intersectTable]=BuildIntersectionTable(testLoop,targLoop,typeLoop)
    
    intersectTable=zeros([numel(testLoop),numel(targLoop)]);
    for ii=1:numel(targLoop)
        for jj=1:numel(testLoop)
            [intersectTable(jj,ii)]=TestInternalOrIntersect(...
                targLoop(ii).coord,testLoop(jj).(typeLoop));
        end
    end
end
%% Test Area error between topo profiles

function [indepLoop]=AreaErrorTopo(testLoop,targLoop)
    % Takes in two sets of loops internally independant and calculates the
    % profiles not filled by loops of both structures
    % This function is recursive and a possible memory hog
    
    % Remove Inner loops from intersect list, their area will simply be
    % substracted at the end.
    [intersectTable]=BuildIntersectionTable(testLoop,targLoop,'coord');
    
    indepTarg=find(all(intersectTable~=-1,1));
    indepTest=find(all(intersectTable~=-1,2));
    indepLoop=repmat(struct('coord',[],'out',[]),[1 numel(indepTarg)+numel(indepTest)]);
    [indepLoop(:).coord]=deal(testLoop(indepTest).coord,targLoop(indepTarg).coord);
    
    for ii=1:numel(indepTest)
        indepLoop(ii).out=(-any(intersectTable(indepTest(ii),:)==2)+0.5)*2;
    end
    for ii=(numel(indepTest)+1):(numel(indepTest)+numel(indepTarg))
        indepLoop(ii).out=(-any(intersectTable(:,indepTarg(ii-numel(indepTest)))==1)+0.5)*2;
    end
    
    % Build groups of tests and targs that intersect.
    
    [testGroups]=BuildGroupLists(intersectTable);
    
    % Intersect each groups
    
    for ii=1:size(testGroups,1)
        targCurr=targLoop(testGroups{ii,1});
        testCurr=testLoop(testGroups{ii,2});
        
        [newloop]=SeparateProfilesArea(testCurr(1).coord,targCurr(1).coord);
        [testCurr]=AreaErrorTopo(newloop,testCurr(2:end));
        [testCurr]=AreaErrorTopo(testCurr,targCurr(2:end));
        
        indepLoop=[indepLoop,testCurr];
    end
    
end

function [testGroups]=BuildGroupLists(intersectTable)
    % Returns groups of mutually intersecting loops
    % Each group is represented by two cells
    % first cell gives the columns in a group
    % second cell gives the rows in a group
    
    actCol=find(any(intersectTable==-1,1));
    actRow=find(any(intersectTable==-1,2))';
    
    kk=0;
    interCols=[];
    testGroups=cell([0 2]);
    while ~isempty(actCol) || ~isempty(actRow)
        
        if isempty(interCols) && (~isempty(actCol) || ~isempty(actRow))
            kk=kk+1;
            interCols=1;
            testGroups(kk,1:2)={[],[]};
        end
        
        
        testGroups{kk,1}=[testGroups{kk,1},actCol(interCols)];
        [tempRows,~]=find(intersectTable(:,actCol(interCols))==-1);
        tempRows=unique(tempRows);
        actCol(interCols)=[];
        interRows=FindObjNum([],tempRows,actRow);
        %tempRows=tempRows(interRows~=0);
        interRows=interRows(interRows~=0);
        
        testGroups{kk,2}=[testGroups{kk,2},actRow(interRows)];
        
        
        [~,tempCols]=find(intersectTable(actRow(interRows),:)==-1);
        tempCols=unique(tempCols);
        actRow(interRows)=[];
        interCols=FindObjNum([],tempCols,actCol);
        %tempCols=tempCols(interCols~=0);
        interCols=interCols(interCols~=0);
    end
    
    %     for ii=1:size(testGroups,1)
    %         testGroups{ii,1}=sort(testGroups{ii,1});
    %         testGroups{ii,2}=sort(testGroups{ii,2});
    %     end
end
%%

function [analysisLoop,targPrep]=PrepareLoopCoord(loop,profileComp,typeLoop)
    % prepares the loop into the right coordinates
    
    nLoop=numel(loop);
    [analysisLoop(1:nLoop).coord]=deal(loop(:).(typeLoop));
    %analysisCoord=vertcat(loop(:).snaxel.coord);
    
    for ii=1:(numel(analysisLoop))
        analysisLoop(ii).nPts=size(analysisLoop(ii).coord,1);
    end
    
    totPts=sum([analysisLoop(:).nPts]);
    
    
    switch profileComp
        case 'distance'
            error('Invalid Profile Comparison method for topology')
        case 'area'
            x=[linspace(1,0,round(totPts/2)),...
                linspace(0,1,round(totPts/2))];
            x=(0.5-cos(x*pi)/2);
            upperLower=ones(size(x));
            upperLower((1:numel(upperLower)<=round(totPts/2)))=-1;
            targPrep=[x',upperLower'];
        case 'distancenorm'
            x=[linspace(1,0,round(totPts/2)),...
                linspace(0,1,round(totPts/2))];
            x=(0.5-cos(x*pi)/2);
            upperLower=ones(size(x));
            upperLower((1:numel(upperLower)<=round(totPts/2)))=-1;
            targPrep=[x',upperLower'];
    end
    
    
    
end


%% Profile extraction

function [nacaCoord]=GenerateLibCoord(x,uplow,airfoilstr)
    % Generates y coordinate for a NACA airfoil.GenerateLibCoord
    % X are the desired x coordinates, uplow is a vector the same size
    
    % xTop and xBot need to be normalised
    
    teps=5.48e-04/2/0.8; % true @ corner=1e-5
    
    
    [airfoilDat,airfoil]=ReadAirfoilData(airfoilstr,'');
    
    xMax=airfoil.func.xMax; %max(xPos);
    xMin=airfoil.func.xMin; %min(xPos);
    x(x>xMax)=xMax;
    x(x<xMin)=xMin;
    
    if any(abs(uplow)~=1)
        error('Vector uplow indicating uppper or lower surface is not well formed')
    end
    y=zeros(size(x));
    y(uplow==1)=airfoil.func.upper(x(uplow==1))+(x(uplow==1)-xMin)/(xMax-xMin)*teps;
    y(uplow==-1)=airfoil.func.lower(x(uplow==-1))-(x(uplow==-1)-xMin)/(xMax-xMin)*teps;
    
    
    if size(x,1)==1
        nacaCoord=[x;y]';
    else
        nacaCoord=[x,y];
    end
    
    
end

function [nacaCoord]=GenerateNacaCoordTOPO(x,uplow,nacaStr)
    % Generates y coordinate for a NACA airfoil.
    % X are the desired x coordinates, uplow is a vector the same size
    
    % xTop and xBot need to be normalised
    a4_open=0.1015;
    a4_closed=0.1036;
    
    naca45t=@(x,t,c,xMin,a4,teps)  5*t*c*(0.2969*sqrt((x-xMin)/c)-0.1260*((x-xMin)/c)...
        -0.3516*((x-xMin)/c).^2+0.2843*((x-xMin)/c).^3-a4*((x-xMin)/c).^4)+((x-xMin)/c)*teps;
    naca4c=@(x,m,p,c,xMin) [m/p^2*(2*p*((x((x-xMin)<(p*c))-xMin)/c)-((x((x-xMin)<(p*c))-xMin)/c).^2),...
        m/(1-p)^2*((1-2*p)+2*p*((x((x-xMin)>=(p*c))-xMin)/c)-((x((x-xMin)>=(p*c))-xMin)/c).^2)];
    teps=5.48e-04/2/0.8; % true @ corner=1e-5
    
    x(x>1)=1;
    x(x<0)=0;
    
    
    
    if numel(nacaStr)==4
%         [ctc,pct,tmax,refFlag]=ReadNacaString(nacaStr);
%         
%         tDist=naca45t(x',tmax,1,0,a4_closed,teps)';
%         [x2,sortOrd]=sort(x);
%         [~,sortOrd2]=sort(sortOrd);
%         cDist=naca4c(x2',ctc,pct,1,0)';
%         cDist=cDist(sortOrd2);
%         if any(abs(uplow)~=1)
%             warning('Vector uplow indicating uppper or lower surface is not well formed')
%         end
%         y=cDist+uplow.*tDist;
%         if size(x,1)==1
%             nacaCoord=[x;y]';
%         else
%             nacaCoord=[x,y];
%         end
        nPtsPloop=size(x,1)/1;
        x=[linspace(0,1,round(nPtsPloop/2))];
        x=(0.5-cos(x*pi)/2);
        [nacaCoord(1).coord,TEPos]=GenerateNACACoordOrient...
                (x,nacaStr,'10','0',...
                [0 0],naca45t,naca4c,a4_closed,teps);
    elseif numel(nacaStr)==5
        [ctc,pct,tmax,refFlag]=ReadNacaString(nacaStr);
        error('Five digits not implemented - need for tabulated data')
    elseif numel(nacaStr)>=8
        % Nine digit and more if for multi element airfoil
        % uses separator ";" between foils
        % uses _ separator between parameters
        % first digit is the number of airfoils
        % then each airfoil is placed as follows:
        % [NACA number l (in 1/10ths of chord) alpha in degrees(LE) xLE (from prev TE)
        % yLE (from prev TE)
        
        
        nacaCell=regexp(nacaStr,';','split');
        nacaCell=regexp(nacaCell,'_','split');
        
        nPtsPloop=size(x,1)/str2double(nacaCell{1}{1});
        
        x=[linspace(0,1,round(nPtsPloop/2))];
        x=(0.5-cos(x*pi)/2);
        TEPos=[0 0];
        for ii=1:str2double(nacaCell{1}{1})
            refPos=TEPos+str2num([nacaCell{ii+1}{4},' ',nacaCell{ii+1}{5}]);
            [nacaCoord(ii).coord,TEPos]=GenerateNACACoordOrient...
                (x,nacaCell{ii+1}{1},nacaCell{ii+1}{2},nacaCell{ii+1}{3},...
                refPos,naca45t,naca4c,a4_closed,teps);
        end
        
        if numel(nacaCell{1})==1
            condNorm='lr';
        else
            condNorm=nacaCell{1}{2};
        end
        
        
        switch condNorm
            case 'lr' % rotate and length normalisation
                normL=sqrt(sum(TEPos.^2));
                rot=-asin(TEPos(2)/normL);
            case 'l'
                normL=sqrt(sum(TEPos.^2));
                rot=0;
            case 'r'
                normL=1;
                rot=-asin(TEPos(2)/normL);
            case 'x'
                normL=abs(TEPos(1));
                rot=0;
        end
        
        
        for ii=1:str2double(nacaCell{1}{1})
            
            nacaCoord(ii).coord=nacaCoord(ii).coord/normL;
            rotMat=[cos(rot),-sin(rot);sin(rot),cos(rot)];
            nacaCoord(ii).coord=(rotMat*(nacaCoord(ii).coord)')';
        end
        
%         plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
%         figure, hold on
%         for ii=1:str2double(nacaCell{1}{1})
%             plotPoints(nacaCoord(ii).coord)
%         end
        
    else
        error('Unrecognised naca String length')
        
    end
    
end

function [coord,TEPos]=GenerateNACACoordOrient(x,naca4Num,l,rot,refPos,tFunc,cFunc,a4,teps)
    
    [ctc,pct,tmax,refFlag]=ReadNacaString(naca4Num);
    
    tDist=tFunc(x,tmax,1,0,a4,teps);
    cDist=cFunc(x,ctc,pct,1,0)';
    
    unitCoord=RemoveIdenticalConsecutivePoints([[x',(tDist'+cDist)];flip([x',(-tDist'+cDist)])]);
    
    unitCoord=unitCoord*str2double(l)/10;
    [~,iTE]=max(unitCoord(:,1));
    rot=-str2double(rot)/180*pi;
    rotMat=[cos(rot),-sin(rot);sin(rot),cos(rot)];
    
    rotCoord=(rotMat*(unitCoord)')';
    
    coord=rotCoord+repmat(refPos,[size(rotCoord,1),1]);
    TEPos=coord(iTE,:);
end

function [targCoord]=LoadLoopData(aeroName)
    
    aeroName=MakePathCompliant(aeroName);
    switch aeroName(end-2:end)
        case 'mat'
            [targCoord]=LoadLoopDataFromMat(aeroName);
        case 'dat'
            targCoord=BoundaryInput(aeroName);
        otherwise
            errstruct.identifier='Optimiser:InvTopo:BadInput:UnrecognisedExtension';
            errstruct.message='for "loop" class of profile to match recognised file extensions are .mat and .dat';
            error(errstruct)
    end
end

function [targCoord]=LoadLoopDataFromMat(aeroName)
    inDat=load(aeroName);
    if isfield(inDat,'loop')
        targCoord=inDat.loop;
        if ~isfield(targCoord,'coord')
            errstruct.identifier='Optimiser:InvTopo:BadInput:LoopNoCoord';
            errstruct.message='Input "loop" variable must have a field called "coord"';
            error(errstruct)
        end
    else
        errstruct.identifier='Optimiser:InvTopo:BadInput:NoLoop';
        errstruct.message=char('.mat file "aeroName"  must contain a variable called "loop"',...
            ['aeroName = ',aeroName]);
        error(errstruct)
    end
end
%% Error Matching

function [errorMeasure,modifiedDistance]=CompareProfilesDistance(profileCoord,targCoord)
    % Compares point to point distance
    
    multipliers=ones(size(targCoord(:,1)));
    multipliers(targCoord(:,1)<0.2*max(targCoord(:,1)))=2;
    
    modifiedDistance=sqrt(sum((profileCoord-targCoord).^2,2)).*multipliers;
    
    errorMeasure.sum=sum(modifiedDistance);
    errorMeasure.mean=mean(modifiedDistance);
    errorMeasure.std=std(modifiedDistance);
    errorMeasure.max=max(modifiedDistance);
    errorMeasure.min=min(modifiedDistance);
    
end

function [errorMeasure,areaDistrib]=CompareProfilesArea(profileCoord,targCoord)
    
    [x0,y0,iout,jout] = intersections(profileCoord([1:end,1],1),profileCoord([1:end,1],2),...
        targCoord([1:end,1],1),targCoord([1:end,1],2));
    rmv=isnan(iout) | isnan(jout);
    iout(rmv)=[];
    x0(rmv)=[];
    y0(rmv)=[];
    jout(rmv)=[];
    
    [iout,sortIndex]=sort(iout);
    x0=x0(sortIndex);
    y0=y0(sortIndex);
    jout=jout(sortIndex);
    targArea=abs(CalculatePolyArea(targCoord));
    nX0=numel(x0);
    nP1=size(profileCoord,1);
    nP2=size(targCoord,1);
    
    if nX0>0
        areaErr=zeros(size(x0))';
        areaLength=zeros(size(x0))';
        areaPosx=zeros(size(x0))';
        areaPosXmin=zeros(size(x0))';
        areaPosXmax=zeros(size(x0))';
        for ii=1:nX0
            % Need to make sure both profiles go in the same direction
            iip1=mod(ii,nX0)+1;
            if floor(iout(ii))~=floor(iout(iip1))
                if ceil(iout(ii))>floor(iout(iip1))
                    ind1=[ceil(iout(ii)):nP1,1:floor(iout(iip1))];
                else
                    ind1=ceil(iout(ii)):floor(iout(iip1));
                end
            else
                ind1=[];
            end
            if floor(jout(ii))~=floor(jout(iip1))
                if ceil(jout(ii))>floor(jout(iip1))
                    ind2=flip([ceil(jout(ii)):nP2,1:floor(jout(iip1))]);
                else
                    ind2=flip(ceil(jout(ii)):floor(jout(iip1)));
                end
            else
                ind2=[];
            end
            
            actPts=[[x0(ii),y0(ii)];
                profileCoord(ind1,:);
                [x0(iip1),y0(iip1)];
                targCoord(ind2,:)];
            areaErr(ii)=abs(CalculatePolyArea(actPts));
            areaLength(ii)=max(actPts(:,1))-min(actPts(:,1));
            areaPosx(ii)=(min(actPts(:,1))+max(actPts(:,1)))/2;
            areaPosXmin(ii)=min(actPts(:,1));
            areaPosXmax(ii)=max(actPts(:,1));
        end
        
        errorMeasure.sum=sum(areaErr)/targArea;
        errorMeasure.mean=mean(areaErr)/targArea;
        errorMeasure.std=std(areaErr)/targArea;
        errorMeasure.max=max(areaErr)/targArea;
        errorMeasure.min=min(areaErr)/targArea;
        
        areaDistrib=[min(areaPosXmin),areaPosx,areaPosXmax;
            0,areaErr./areaLength*2,-areaErr./areaLength*2];
        [areaDistrib(1,:),iSortDistrib]=sort(areaDistrib(1,:));
        eqInd=find(abs(areaDistrib(1,1:end-1)-areaDistrib(1,2:end))<1e-10);
        areaDistrib(2,:)=(areaDistrib(2,iSortDistrib));
        areaDistrib(2,eqInd)=areaDistrib(2,eqInd)+areaDistrib(2,eqInd+1);
        areaDistrib(:,eqInd+1)=[];
        areaDistrib(2,:)=cumsum(areaDistrib(2,:));
        areaDistrib(:,find(areaDistrib(2,:)<1e-18))=[];
        %areaDistrib(:,find(areaDistrib(1,1:end-1)<1e-15))=[];
        areaDistrib=areaDistrib';
        %         distribError=areaErr./areaLength;
        %         xPts=[[profileCoord(1,1),reshape(x0,[1,numel(x0)])];
        %             [reshape(x0,[1,numel(x0)]),profileCoord(end,1)]];
        %         xPts=xPts(:);
        %         errPts=[distribError(end),distribError;distribError(end),distribError];
        %         errPts=errPts(:);
        %         errorDist=[xPts',errPts'];
    else
        
        [errorMeasure2]=abs(abs(CalculatePolyArea(profileCoord))-targArea);
        errorMeasure.sum=errorMeasure2/targArea;
        errorMeasure.mean=errorMeasure.sum;
        errorMeasure.std=errorMeasure.sum;
        errorMeasure.max=errorMeasure.sum;
        errorMeasure.min=errorMeasure.sum;
        areaDistrib=[(min(profileCoord(:,1))+min(profileCoord(:,1)))/2,errorMeasure.sum];
    end
    
    %error('Compare Profile through area integration has not been coded yet')
end


function [errorMeasure,areaDistrib,loop]=IndepProfileError(loop,targArea,funArea)
    if nargin<3
        funArea{1}=@(x) x;
        funArea{2}=@(x) x;
        if nargin<2
            targArea=1;
        end
    end
    
    nX0=numel(loop);
    
    areaErr=zeros(size(loop));
    areaLength=zeros(size(loop));
    areaPosx=zeros(size(loop));
    areaPosXmin=zeros(size(loop));
    areaPosXmax=zeros(size(loop));
    for ii=1:nX0
        
        
        actPts=loop(ii).coord;
        loop(ii).normA=funArea{1}(abs(CalculatePolyArea(actPts))/targArea);
        areaErr(ii)=loop(ii).out*loop(ii).normA;
        areaLength(ii)=max(actPts(:,1))-min(actPts(:,1));
        areaPosx(ii)=(min(actPts(:,1))+max(actPts(:,1)))/2;
        areaPosXmin(ii)=min(actPts(:,1));
        areaPosXmax(ii)=max(actPts(:,1));
    end
    
    errorMeasure.sum=funArea{2}(sum(areaErr));
    errorMeasure.mean=funArea{2}(mean(areaErr));
    errorMeasure.std=funArea{2}(std(areaErr));
    errorMeasure.max=funArea{2}(max(areaErr));
    errorMeasure.min=funArea{2}(min(areaErr));
    
    areaDistrib=[min(areaPosXmin),areaPosx,areaPosXmax;
        0,areaErr./areaLength*2,-areaErr./areaLength*2];
    [areaDistrib(1,:),iSortDistrib]=sort(areaDistrib(1,:));
    eqInd=find(abs(areaDistrib(1,1:end-1)-areaDistrib(1,2:end))<1e-10);
    areaDistrib(2,:)=(areaDistrib(2,iSortDistrib));
    areaDistrib(2,eqInd)=areaDistrib(2,eqInd)+areaDistrib(2,eqInd+1);
    areaDistrib(:,eqInd+1)=[];
    areaDistrib(2,:)=cumsum(areaDistrib(2,:));
    %areaDistrib(:,find(areaDistrib(2,:)<1e-18))=[];
    %areaDistrib(:,find(areaDistrib(1,1:end-1)<1e-15))=[];
    areaDistrib=areaDistrib';
    %         distribError=areaErr./areaLength;
    %         xPts=[[profileCoord(1,1),reshape(x0,[1,numel(x0)])];
    %             [reshape(x0,[1,numel(x0)]),profileCoord(end,1)]];
    %         xPts=xPts(:);
    %         errPts=[distribError(end),distribError;distribError(end),distribError];
    %         errPts=errPts(:);
    %         errorDist=[xPts',errPts'];
    
end

function [newloop]=SeparateProfilesArea(profileCoord,targCoord)
    % compares two loops and returns the area not filled by both
    if ~CCWLoop(profileCoord)
        profileCoord=flip(profileCoord);
    end
    if ~CCWLoop(targCoord)
        targCoord=flip(targCoord);
    end
    
    [x0,y0,iout,jout] = intersections(profileCoord([1:end,1],1),profileCoord([1:end,1],2),...
        targCoord([1:end,1],1),targCoord([1:end,1],2));
    rmv=isnan(iout) | isnan(jout);
    iout(rmv)=[];
    x0(rmv)=[];
    y0(rmv)=[];
    jout(rmv)=[];
    
    [iout,sortIndex]=sort(iout);
    x0=x0(sortIndex);
    y0=y0(sortIndex);
    jout=jout(sortIndex);
    
    nX0=numel(x0);
    nP1=size(profileCoord,1);
    nP2=size(targCoord,1);
    
    if nX0>0
        newloop=repmat(struct('coord',zeros(0,2)),[1,nX0]);
        nXplore=1:nX0;
        pLoop=0;
        while ~isempty(nXplore)
            % Need to make sure both profiles go in the same direction
            
            % THis code is insufficient it does not check for that no other
            % intersection lies on ind2 (cannot lie on ind1 as it is sorted)
            kk=1;
            iiStart=nXplore(kk);
            condLoop=true;
             pLoop=pLoop+1;
            while condLoop
                ii=nXplore(kk);
                nXplore(kk)=[];
                iip1=mod(ii,nX0)+1;
                if floor(iout(ii))~=floor(iout(iip1))
                    if ceil(iout(ii))>floor(iout(iip1))
                        ind1=[ceil(iout(ii)):nP1,1:floor(iout(iip1))];
                    else
                        ind1=ceil(iout(ii)):floor(iout(iip1));
                    end
                else
                    ind1=[];
                end
                iij=ii;
                if floor(jout(ii))~=floor(jout(iip1))
                    if ceil(jout(ii))>floor(jout(iip1))
                        
                        isInTheWay1=(0<jout & jout<jout(iip1));
                        isInTheWay2=(jout>jout(ii) & jout<nP2);
                        if any(isInTheWay1)
                            subIsInTheWay=find(isInTheWay1);
                            [~,posInWay]=max(jout(subIsInTheWay));
                            iij=subIsInTheWay(posInWay);
                            ind2=flip(ceil(jout(iij)):floor(jout(iip1)));
                        elseif any(isInTheWay2)
                            subIsInTheWay=find(isInTheWay2);
                            [~,posInWay]=max(jout(subIsInTheWay));
                            iij=subIsInTheWay(posInWay);
                            ind2=flip([ceil(jout(iij)):nP2,1:floor(jout(iip1))]);
                        else
                            ind2=flip([ceil(jout(iij)):nP2,1:floor(jout(iip1))]); % iip1 does not change
                        end
                    else
                        
                        isInTheWay=jout>jout(ii) & jout<jout(iip1);
                        if any(isInTheWay)
                            subIsInTheWay=find(isInTheWay);
                            [~,posInWay]=max(jout(subIsInTheWay));
                            iij=subIsInTheWay(posInWay);
                        end
                        ind2=flip(ceil(jout(iij)):floor(jout(iip1)));
                        
                    end
                else
                    ind2=[];
                end
                
                
                
                newloop(pLoop).coord=[newloop(pLoop).coord;
                    [x0(ii),y0(ii)];
                    profileCoord(ind1,:);
                    [x0(iip1),y0(iip1)];
                    targCoord(ind2,:)];
                kk=find(nXplore==iij);
                condLoop=iij~=iiStart && ~isempty(kk);
            end
           
            
        end
         
        newloop=newloop(1:pLoop);
        keepLoop = true(size(newloop));
        
        for ii=1:pLoop
            keepLoop(ii) = size(newloop(ii).coord,1)>2;
        end
        newloop=newloop(keepLoop);
        
        
        
    else
        newloop=repmat(struct('coord',zeros(0,2)),[1,1]);
%         newloop.coord = profileCoord;
        error('No Intersection, should not happen')
        % can happen as the table intersection is not the same as the
        % actual intersection
    end
    
    %error('Compare Profile through area integration has not been coded yet')
end

function [errorMeasure,areaDistrib]=CompareProfilesDistance2(profileCoord,targCoord)
    
    
    normVec=@(v) v./repmat(sqrt(sum(v.^2,2)),[1,size(v,2)]);
    pDistDir=@(m,p,v) sign((m-repmat(p,[size(m,1),1]))*v').*sqrt(sum((m-repmat(p,[size(m,1),1])).^2,2));
    
    
    
    maxDist=[min([profileCoord;targCoord]);max([profileCoord;targCoord])];
    nP1=size(profileCoord,1);
    edgeNormal=([0 1;-1 0]*normVec(profileCoord([2:end,1],:)-profileCoord)')';
    [~,vDist]=knnsearch(targCoord,profileCoord);
    for ii=1:nP1
        
        vertNormal=sum(edgeNormal(mod((ii-1:ii)-1,nP1)+1,:));
        
        scales=(maxDist-profileCoord([ii,ii],:))./repmat(vertNormal,[2,1]);
        scales=scales(:);
        scales=[min(scales(scales>0)),max(scales(scales<0))];
        
        testCoord=repmat(profileCoord(ii,:),[2,1])+(scales'*vertNormal);
        
        [x0,y0,iout,jout] = intersections(testCoord(:,1),testCoord(:,2),...
            targCoord([1:end,1],1),targCoord([1:end,1],2));
        rmv=isnan(iout) | isnan(jout);
        
        x0(rmv)=[];
        y0(rmv)=[];
        
        
        if numel(x0)>0
            vertDists=pDistDir([x0,y0],profileCoord(ii,:),vertNormal);
            
            vDist(ii)=min([vDist(ii);abs(vertDists)]);
            
        end
    end
    
    
    errorMeasure.sum=sum(vDist);
    errorMeasure.mean=mean(vDist);
    errorMeasure.std=std(vDist);
    errorMeasure.max=max(vDist);
    errorMeasure.min=min(vDist);
    areaDistrib=[profileCoord(:,1),vDist];
    
    
    %error('Compare Profile through area integration has not been coded yet')
end














