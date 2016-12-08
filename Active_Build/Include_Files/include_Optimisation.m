function [] = include_Optimisation()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%%
function [constrVal]=NacaOuterLimit0012(gridrefined,paramoptim)
    
    [constrVal]=NacaOuterLimit4d(gridrefined,paramoptim,'0012');
    %{
    naca4t=@(x,t,c,xMin)  5*t*c*(0.2969*sqrt((x-xMin)/c)-0.1260*((x-xMin)/c)...
        -0.3516*((x-xMin)/c).^2+0.2843*((x-xMin)/c).^3-0.1036*((x-xMin)/c).^4);
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    
    t=0.12;
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
    actVerts=[cellCentredGrid(isActive).vertex];
    coord=vertcat(actVerts(:).coord);
    
    coord(:,2)=coord(:,2)*axisRatio;%;
    
    xPos=RemoveIdenticalEntries(coord(:,1));
    xMax=max(xPos);
    xMin=min(xPos);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        x=linspace(posMin(1),posMax(1),200);
        tDistrib=naca4t(x,t,(xMax-xMin),xMin);
        y=min(max(tDistrib,posMin(2)),posMax(2))-min(max(-tDistrib,posMin(2)),posMax(2));
        %plot(x,y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
    %}
end

function [constrVal]=NacaOuterLimit4d(gridrefined,paramoptim,nacaStr)
    a4_open=0.1015;
    a4_closed=0.1036;
    naca4t=@(x,t,c,xMin,a4,teps)  5*t*c*(0.2969*sqrt((x-xMin)/c)-0.1260*((x-xMin)/c)...
        -0.3516*((x-xMin)/c).^2+0.2843*((x-xMin)/c).^3-a4*((x-xMin)/c).^4)+((x-xMin)/c)*teps;
    
    naca4c=@(x,m,p,c,xMin) [m/p^2*(2*p*((x((x-xMin)<(p*c))-xMin)/c)-((x((x-xMin)<(p*c))-xMin)/c).^2),...
        m/(1-p)^2*((1-2*p)+2*p*((x((x-xMin)>=(p*c))-xMin)/c)-((x((x-xMin)>=(p*c))-xMin)/c).^2)];
%     ((x(x>=(p*c))-xMin)/c)
%     
%     [m*x(x<(cp*c))/p^2.*(2*p-x(x<(p*c))) ;...
%         m*(1-x(x>=(p*c))/c)/(1-p)^2.*(1+x(x>=(p*c))-2*p)];
    
    teps=5.48e-04/2/0.8;
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    [m,p,t,refFlag]=ReadNacaString(nacaStr);
    if m==0 || p==0
        p=0.0001;
        m=0;
    end
    
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
    actVerts=[cellCentredGrid(isActive).vertex];
    coord=vertcat(actVerts(:).coord);
    %figure, hold on
    coord(:,2)=coord(:,2)*axisRatio;%;
    
    xPos=RemoveIdenticalEntries(coord(:,1));
    xMax=1; %max(xPos);
    xMin=0; %min(xPos);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        %x=linspace(posMin(1),posMax(1),200);
        x=min(max(linspace(posMin(1),posMax(1),200),xMin),xMax);
        tDistrib=naca4t(x,t,(xMax-xMin),xMin,a4_closed,teps);
        cDistrib=naca4c(x,m,p,(xMax-xMin),xMin);
        y=min(max(cDistrib+tDistrib,posMin(2)),posMax(2))-min(max(cDistrib-tDistrib,posMin(2)),posMax(2));
        %plot(x,cDistrib+y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
end

function [airfoilDat,airfoil]=ReadAirfoilData(airfoilstr,airfoilDir)
    
    if nargin==1 || isempty(airfoilDir)
        airfoilDir='..\AeroLib\Smoothed\';
    end
    
    fid=fopen([airfoilDir,airfoilstr,'.dat'],'r');
    if fid==-1
        error('Airfoil Data could not be read. %s was not recognised as a valid file',[airfoilDir,airfoilstr,'.dat']);
    end
    for ii=1:3, str=fgetl(fid);end
    
    [airfoilDat,check]=fscanf(fid,'%f %f \n',[2,inf]);
    airfoilDat=airfoilDat';
    fclose(fid);
    
    [isCCW]=CCWLoop(airfoilDat);
    if ~isCCW
        airfoilDat=flip(airfoilDat,1);
    end
    iTE=find(airfoilDat(:,1)==max(airfoilDat(:,1)));
    if numel(iTE)==1
        airfoilDat=airfoilDat([iTE:end,1:iTE],:); % repeats the trailing edge at start and end.
    else
        [~,teT]=max(airfoilDat(iTE,2));
        [~,teB]=min(airfoilDat(iTE,2));
        if iTE(teT)==1
            airfoilDat=airfoilDat([iTE(teT):iTE(teB)],:);
        else
            airfoilDat=airfoilDat([iTE(teT):end,1:iTE(teB)],:);
        end
    end
    %analysisCoord(:,2)=analysisCoord(:,2)-analysisCoord(1,2); % 0 the y at the trailing edge.
    % TE is start and end of coord list, need LE
    %[dLE,iLE]=max(sum((analysisCoord-ones([size(analysisCoord,1),1])*analysisCoord(1,:)).^2,2));
    [~,iLE]=min(airfoilDat(:,1));
    if numel(iLE)==1
        airfoilDat=airfoilDat([1:iLE,iLE:end],:); % repeats the trailing edge at start and end.
        leT=1;
        leB=2;
        iLE=[iLE,iLE+1];
    else
        [~,leT]=max(airfoilDat(iLE,2));
        [~,leB]=min(airfoilDat(iLE,2));
        airfoilDat=airfoilDat([1:iLE(leT),iLE(leB):end],:);
    end
    
    
    upperLower=ones(size(airfoilDat(:,1)));
    teps=ones(size(airfoilDat(:,1)))*airfoilDat(iTE(teT),2);
    upperLower(iLE(leB):end)=-1;
    teps(iLE(leB):end)=airfoilDat(iTE(teB),2);
    airfoil.func.xMin=min(airfoilDat(:,1));
    airfoil.func.xMax=max(airfoilDat(:,1));
    airfoilDat(:,2)=airfoilDat(:,2)-teps.*(airfoilDat(:,1)-airfoil.func.xMin)...
        /(airfoil.func.xMax-airfoil.func.xMin);
    
    airfoil.coord=airfoilDat;
    airfoil.upperlower=upperLower;
    airfoil.func.upper=griddedInterpolant(flip(airfoilDat(1:iLE(leT),1)),...
        flip(airfoilDat(1:iLE(leT),2)),'linear','none');
    airfoil.func.lower=griddedInterpolant(airfoilDat(iLE(leB):end,1),...
        airfoilDat(iLE(leB):end,2),'linear','none');
end

function [constrVal]=OuterLimitInverse(gridrefined,paramoptim,airfoilstr)
    

   
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
     teps=5.48e-04/2/0.8;
     
     
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    [airfoilDat,airfoil]=ReadAirfoilData(airfoilstr,'');
    
    
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
%     actVerts=[cellCentredGrid(isActive).vertex];
%     coord=vertcat(actVerts(:).coord);
%     figure, hold on
%     coord(:,2)=coord(:,2)*axisRatio;%;
    
    xMax=airfoil.func.xMax; %max(xPos);
    xMin=airfoil.func.xMin; %min(xPos);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        %x=linspace(posMin(1),posMax(1),200);
        x=min(max(linspace(posMin(1),posMax(1),200),xMin),xMax);
        uppersurf=airfoil.func.upper(x)+(x-xMin)/(xMax-xMin)*teps;
        lowersurf=airfoil.func.lower(x)-(x-xMin)/(xMax-xMin)*teps;
        y=min(max(uppersurf,posMin(2)),posMax(2))-min(max(lowersurf,posMin(2)),posMax(2));
        %plot(x,cDistrib+y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
end

function [constrVal]=MatchVoltoShape(gridrefined,paramoptim,shapePath)
    
    
%     ((x(x>=(p*c))-xMin)/c)
%     
%     [m*x(x<(cp*c))/p^2.*(2*p-x(x<(p*c))) ;...
%         m*(1-x(x>=(p*c))/c)/(1-p)^2.*(1+x(x>=(p*c))-2*p)];
    
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    [uppersurf,lowersurf]=ReadShapeIn(shapePath);
   
    
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
    actVerts=[cellCentredGrid(isActive).vertex];
    coord=vertcat(actVerts(:).coord);
    %figure, hold on
    coord(:,2)=coord(:,2)*axisRatio;%;
    
    xPos=RemoveIdenticalEntries(coord(:,1));
    xMax=max(xPos);
    xMin=min(xPos);
    xDist=(xMax-xMin);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        x=linspace(posMin(1),posMax(1),200);
        
        yup=interp1(uppersurf(:,1)*xDist+xMin,uppersurf(:,2)*xDist,x);
        ylow=interp1(lowersurf(:,1)*xDist+xMin,lowersurf(:,2)*xDist,x);
        y=min(max(yup,posMin(2)),posMax(2))-min(max(ylow,posMin(2)),posMax(2));
        %plot(x,cDistrib+y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
end

function [uppersurf,lowersurf]=ReadShapeIn(shapepath)
    
    filename=dir(shapepath);
    filename=filename(1).name;
    extPos=regexp(filename,'\.');
    ext=filename(extPos+1:end);
    
    switch ext
        case 'mat'
            load(shapepath)
            if ~exist('uppersurf','var');error('Data file loaded did not have the upper surface');end
            if ~exist('lowersurf','var');error('Data file loaded did not have the lower surface');end
            
        case 'dat'
            error('not coded')
            
        otherwise
            error('Unknown type')
            
    end
        
end

function [ctc,pct,ttc,refFlag]=ReadNacaString(nacaStr)
    refFlag=0;
    if numel(nacaStr)==4
        ctc=str2num(nacaStr(1))/100;
        pct=str2num(nacaStr(2))/10;
        ttc=str2num(nacaStr(3:4))/100;
        
        
    elseif numel(nacaStr)==5
        ctc=str2num(nacaStr(1))/100;
        pct=str2num(nacaStr(2))/10;
        refFlag=str2num(nacaStr(3));
        ttc=str2num(nacaStr(4:5))/100;
        %error('Five digits not implemented - need for tabulated data')
    else
        error('Invalid string length')
    end
end

function [cornersCoord]=IdentifyCorners(coord)
    
    minXPos=find(coord(:,1)==min(coord(:,1)));
    maxXPos=find(coord(:,1)==max(coord(:,1)));
    minYPos=find(coord(:,2)==min(coord(:,2)));
    maxYPos=find(coord(:,2)==max(coord(:,2)));
    
    [~,loleCorn]=min(coord(minXPos,2));
    loleCorn=minXPos(loleCorn);
    [~,hileCorn]=max(coord(minXPos,2));
    hileCorn=minXPos(hileCorn);
    [~,loriCorn]=min(coord(maxXPos,2));
    loriCorn=maxXPos(loriCorn);
    [~,hiriCorn]=max(coord(maxXPos,2));
    hiriCorn=maxXPos(hiriCorn);
    
    
    [~,leloCorn]=min(coord(minYPos,1));
    leloCorn=minYPos(leloCorn);
    [~,riloCorn]=max(coord(minYPos,1));
    riloCorn=minYPos(riloCorn);
    [~,lehiCorn]=min(coord(maxYPos,1));
    lehiCorn=maxYPos(lehiCorn);
    [~,rihiCorn]=max(coord(maxYPos,1));
    rihiCorn=maxYPos(rihiCorn);
    
    tests=(loleCorn==leloCorn) && (loriCorn==riloCorn) &&(hileCorn==lehiCorn) && (hiriCorn==rihiCorn);
    if tests==false
        warning('Unexpected Corners returning only one set')
    end
    cornersCoord=coord([loleCorn,loriCorn,hileCorn,hiriCorn],:);
end

function [desVarList]=ExtractActiveVariable(nFill,notDesInd,inactiveVar)
    
    desVarList=1:nFill;
    desVarList([notDesInd,inactiveVar])=[];
    
end

function [inactiveVar]=SelectInactiveVariables(newFill,varActive,derivtenscalc,cutoff)
    
    switch varActive
        case 'all'
            inactiveVar=[];
        case 'border'
            [inactiveVar]=InactiveVariables_border(newFill);
            
        case 'wideborder'
            [inactiveVar,activeVar]=InactiveVariables_wideborder(newFill,derivtenscalc,cutoff);
            
        case 'snaksensiv'
            
            [inactiveVar]=InactiveVariables_border(newFill);
        otherwise
            error('unrecognised variable activation criterion')
    end
    
    
end

function [inactiveVar,activeVar]=InactiveVariables_border(newFill)
    
    is0=newFill==0;
    is1=newFill==1;
    
    inactiveVar=find(is0);
    inactiveVar=[inactiveVar,find(is1)];
    
    isBord=~is0 & ~is1;
    activeVar=find(isBord);
    
    
end

function [inactiveVar,activeVar]=InactiveVariables_wideborder(newFill,derivtenscalc,cutoff)
    
    [inactiveVar,activeVar]=InactiveVariables_border(newFill);
    
    indexDeriv=[derivtenscalc(:).index];
    actVarSub=FindObjNum([],activeVar,indexDeriv);
    actFill=newFill(activeVar);
    cutActVar=actFill>cutoff;
    activeVar=unique([activeVar,[derivtenscalc(actVarSub(cutActVar)).neighbours]]);
    inactiveVar=1:length(newFill);
    inactiveVar(activeVar)=[];
end

function [isGradient]=CheckIfGradient(optimMethod)
    
    switch optimMethod
        
        case 'DE'
            isGradient=false;
        case 'DEtan'
            isGradient=false;
        case 'conjgrad'
            isGradient=true;
        case 'conjgradls'
            isGradient=true;
        otherwise
            isGradient=false;
            warning('Optimisation method is not known as gradient based or otherwise, no gradient is assumed')
            
    end
end

function [isAnalytical]=CheckIfAnalytical(objectiveName)
    
    switch objectiveName
        
        case 'Rosenbrock'
            isAnalytical=true;
        otherwise
            isAnalytical=false;
            warning('Objective is not specified assumed non analytical function')
            
    end
end

function [isSnakeSensitivity]=CheckSnakeSensitivityAlgorithm(paramoptim)
    
    varExtract={'varActive','optimMethod'};
    [varActive,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    isSnakeSensitivity=false;
    [isGradient]=CheckIfGradient(optimMethod);
    
    if isGradient
        isSnakeSensitivity=strcmp(varActive,'snaksensiv');
    end
end

function [newRootFill]=OverflowHandling(paramoptim,newRootFill)
    % Function which handles steps overflowing the fill constraint
    
    varExtract={'desVarRange','varOverflow'};
    [desVarRange,varOverflow]=ExtractVariables(varExtract,paramoptim);
    
    switch varOverflow
        case 'truncate'
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
            
        case 'spill'
            newFill=zeros(size(newRootFill));
            for ii=1:size(newRootFill,1)
                [newFill(ii,:)]=SpillOverflow(paramoptim,newRootFill(ii,:));
            end
            newRootFill=newFill;
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
        otherwise
            error('No valid design variable overflow mechanism')
            
    end
    
    
    
end

function [newFill]=SpillOverflow(paramoptim,newRootFill)
    
    varExtract={'desVarRange','desvarconnec'};
    [desVarRange,desvarconnec]=ExtractVariables(varExtract,paramoptim);
    
    normRootFill=(newRootFill-min(desVarRange))/(max(desVarRange)-min(desVarRange));
    
    overVar=find(normRootFill>1);
    underVar=find(normRootFill<0);
    exitFlag=false;
    kk=0;
    while (~isempty(overVar) || ~isempty(underVar)) && ~exitFlag
        
        
        [normRootFill]=SpillOverflowVarHandling(normRootFill,desvarconnec,overVar);
        
        normRootFill=1-normRootFill;
        [normRootFill]=SpillOverflowVarHandling(normRootFill,desvarconnec,underVar);
        normRootFill=1-normRootFill;
        
        
        
        overVar=find(normRootFill>1);
        underVar=find(normRootFill<0);
        meanFill=mean(normRootFill);
        stdFill=std(normRootFill);
        kk=kk+1;
        if ((meanFill>=1 || meanFill<=0)  && (abs(stdFill)<1e-5)) || kk>100
            exitFlag=true;
        end
        
    end
    
    newFill=normRootFill*(max(desVarRange)-min(desVarRange))+min(desVarRange);
end

function [newRootFill]=SpillOverflowVarHandling(newRootFill,desvarconnec,flowVar)
    
    desVarIndList=[desvarconnec(:).index];
    
    overFlowMat=zeros([length(flowVar),length(newRootFill)]);
    
    for ii=1:length(flowVar)
        
        currSub=FindObjNum([],flowVar(ii),desVarIndList);
        
        neighInd=desvarconnec(currSub).neighbours;
        cornInd=desvarconnec(currSub).corners;
        
        neighSub=FindObjNum([],neighInd,desVarIndList);
        cornSub=FindObjNum([],cornInd,desVarIndList);
        
        currVol=newRootFill(currSub);
        neighVol=newRootFill(neighSub);
        cornVol=newRootFill(cornSub);
        
        neighEmpt=neighSub(neighVol<=0);
        cornEmpt=cornSub(cornVol<=0);
        neighGreyCell={[]};
        for jj=1:length(cornEmpt)
            neighGreyCell{jj}=FindObjNum([],desvarconnec(cornEmpt(jj)).neighbours,neighInd)';
        end
        neighAll=[neighGreyCell{:}];
        neighSubAll=neighSub(RemoveIdenticalEntries(neighAll(neighAll~=0)));
        neighSubAll=RemoveIdenticalEntries([neighEmpt;neighSubAll]);
        neighSubAll=neighSubAll(newRootFill(neighSubAll)<1);
        
        nCorn=numel(cornEmpt);
        nNeigh=sum(1-newRootFill(neighSubAll));
        isNormProb=true;
        
        if nCorn>0
            baseRate=(-nNeigh+sqrt(nNeigh^2+4*nCorn))/(2*nCorn);
            
        elseif nNeigh>0
            baseRate=1/nNeigh;
            
        elseif numel(neighSub(newRootFill(neighSub)<1))>0 ...
                && numel(cornSub(newRootFill(cornSub)<1))>0
            
            neighSubAll=neighSub(newRootFill(neighSub)<1);
            cornEmpt=cornSub(newRootFill(cornSub)<1);
            nCorn=numel(cornEmpt);
            nNeigh=sum(1-newRootFill(neighSubAll));
            
            baseRate=(-nNeigh+sqrt(nNeigh^2+4*nCorn))/(2*nCorn);
            
        elseif numel(neighSub(newRootFill(neighSub)<1))>0
            neighSubAll=neighSub(newRootFill(neighSub)<1);
            nNeigh=sum(1-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
        elseif numel(cornSub(newRootFill(cornSub)<1))>0
            neighSubAll=cornSub(newRootFill(cornSub)<1);
            nNeigh=sum(1-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
        elseif numel(neighSub(newRootFill(neighSub)<currVol))>0 || ...
                numel(cornSub(newRootFill(cornSub)<currVol))>0
            neighSubAll=[neighSub(newRootFill(neighSub)<currVol);cornSub(newRootFill(cornSub)<currVol)];
            nNeigh=sum(currVol-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
            isNormProb=false;
        else
            neighSubAll=[];
            nNeigh=sum(-newRootFill(neighSubAll));
            
            baseRate=0;
            
            isNormProb=false;
        end
        
        if isNormProb
            overFlowVol=zeros(size(newRootFill));
            overFlowVol(neighSubAll)=(1-newRootFill(neighSubAll))*baseRate;
            overFlowVol(cornEmpt)=baseRate^2;
            overFlowVol(currSub)=-1;
            overFlowMat(ii,:)=overFlowVol*(currVol-1);
            %newRootFill=newRootFill+overFlowVol*(currVol-1);
        else
            overFlowVol=zeros(size(newRootFill));
            overFlowVol(neighSubAll)=baseRate*(currVol-newRootFill(neighSubAll));
            overFlowVol(cornEmpt)=baseRate^2;
            overFlowVol(currSub)=-sum(overFlowVol);
            overFlowMat(ii,:)=overFlowVol*(currVol-max([min(newRootFill(neighSubAll)),1]))*2/3;
            %newRootFill=newRootFill+overFlowVol*(currVol-max([min(newRootFill(neighSubAll)),1]))*2/3;
        end
        
    end
    
    newRootFill=newRootFill+sum(overFlowMat,1);
    
end

function population=ApplySymmetry(paramoptim,population)
    
    varExtract={'symDesVarList'};
    [symDesVarList]=ExtractVariables(varExtract,paramoptim);
    
    for ii=1:length(population)
        population(ii).fill(symDesVarList(2,:))=...
            population(ii).fill(symDesVarList(1,:));
    end
    
    
end

function [res]=TestDeviationPopStruct(population,cellTest,dirOut)
    
    for ii=1:numel(cellTest)
        analysisVal=[];
        for jj=1:numel(population)
            analysisVal=[analysisVal,population(jj).additional.(cellTest{ii})];
        end
        h(ii)=figure('Name',cellTest{ii});
        subplot(2,2,1)
        normplot(analysisVal)
        res.(cellTest{ii}).mean=mean(analysisVal);
        res.(cellTest{ii}).median=median(analysisVal);
        res.(cellTest{ii}).std=std(analysisVal);
        res.(cellTest{ii}).max=max(analysisVal);
        res.(cellTest{ii}).min=min(analysisVal);
        res.(cellTest{ii}).decile=quantile(analysisVal,19);
        res.(cellTest{ii}).ratiosuccess=sum(analysisVal<(8e-4))/numel(analysisVal);
        subplot(2,2,3)
        semilogy(res.(cellTest{ii}).decile)
        grid on
        xlabel('decile')
        ylabel('Val')
        analysisVal=log10(analysisVal);
        subplot(2,2,2)
        normplot(analysisVal)
        res.(cellTest{ii}).logmean=10^mean(analysisVal);
        res.(cellTest{ii}).logmedian=10^median(analysisVal);
        res.(cellTest{ii}).logstd=std(analysisVal);
        res.(cellTest{ii}).logmax=10^max(analysisVal);
        res.(cellTest{ii}).logmin=10^min(analysisVal);
        res.(cellTest{ii}).logdecile=10.^quantile(analysisVal,19);
        subplot(2,2,4)
        plot(quantile(analysisVal,19))
        grid on
        xlabel('decile')
        ylabel('log10(Val)')
    end
    
    if nargin>2
        for ii=1:numel(h)
            figName=[dirOut,filesep,'Figure_',int2str(ii),'_',h(ii).Name,'.fig'];
            hgsave(h(ii),figName);
        end
        fileOut=fopen([dirOut,filesep,'TabResSummary.csv'],'w');
        fieldsRes=fieldnames(res.(cellTest{1}));
        
        fprintf(fileOut,' ');
        for ii=1:numel(fieldsRes)
            fprintf(fileOut,', %s',fieldsRes{ii});
        end
        
        for jj=1:numel(cellTest)
            fprintf(fileOut,'\n');
            fprintf(fileOut,'%s',cellTest{jj});
            for ii=1:numel(fieldsRes)
                fprintf(fileOut,',%s',...
                    num2str(res.(cellTest{jj}).(fieldsRes{ii}),'%20.8e '));
            end
        end
        fclose(fileOut);
    end
end
