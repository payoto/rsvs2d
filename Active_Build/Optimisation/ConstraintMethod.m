%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%            Optimisation Using
%          Parametric Snakes for
%         for Aerodynamic shape
%         parametrisation
%             Constraint Handling
%
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [population,varargout]=ConstraintMethod(entryPoint,paramoptim,population,varargin)
    % Function distributing the optimisation to various methods
    
    varargout{1}=[];
    switch entryPoint
        case 'init'
            varExtract={'initConstr','initVal'};
            [initConstr,initVal]=ExtractVariables(varExtract,paramoptim);
            
            for ii=1:length(initConstr)
                
                [population,varargout]=InitVariableConsCaller(initConstr{ii},initVal{ii},...
                    paramoptim,population,varargin{:});
            end
            if isempty(varargout{1}),varargout{1}=paramoptim;end
            
        case 'DesVar'
            varExtract={'desVarConstr','desVarVal'};
            [desVarConstr,desVarVal]=ExtractVariables(varExtract,paramoptim);
            
            for ii=1:length(desVarConstr)
                
                [population]=DesignVariableConsCaller(desVarConstr{ii},desVarVal{ii},...
                    paramoptim,population,varargin{:});
            end
            
        case 'Res'
            varExtract={'resConstr','resVal'};
            [resConstr,resVal]=ExtractVariables(varExtract,paramoptim);
            
            for ii=1:length(resConstr)
                [population]=ResultVariableConsCaller(resConstr{ii},resVal{ii},...
                    paramoptim,population,varargin{:});
            end
    end
    
end

%% Initialisation cases

function [population,argout]=InitVariableConsCaller(constrName,constrVal,paroptim,population,varargin)
    argout{1}=[];
    switch constrName
        case 'LocalVolFrac_image'
            [paramoptim]=LocalConstraintExtraction(paroptim,constrVal);
            argout{1}=paramoptim;
        case 'Naca0012'
            
        case ' '
            
        otherwise
            error('Design Variable Constraint Not Recognised')
    end
    
    
    
end

function [paramoptim]=LocalConstraintExtraction(paramoptim,constrVal)
    
    varExtract={'cellLevels','desVarConstr','desVarVal'};
    [cellLevels]=ExtractVariables(varExtract(1),paramoptim.parametrisation);
    [desVarConstr,desVarVal]=ExtractVariables(varExtract(2:3),paramoptim);
    
    [desVal,cellLevels]=ImageProcess(constrVal{1},cellLevels);
    desVarConstr{end+1}=['LocalVolFrac_',constrVal{2}];
    desVarVal{end+1}=desVal;
    
    varExtract={'cellLevels','passDomBounds','desVarConstr','desVarVal'};
    [paramoptim.parametrisation]=SetVariables(varExtract(1:2),{cellLevels,...
        MakeCartesianGridBounds(cellLevels)},paramoptim.parametrisation);
    [paramoptim.initparam]=SetVariables(varExtract(1:2),{cellLevels,...
        MakeCartesianGridBounds(cellLevels)},paramoptim.initparam);
    [paramoptim]=SetVariables(varExtract(3:4),{desVarConstr,desVarVal},paramoptim);
end

function [constrVal,cellLevels]=ImageProcess(imPath,cellLevels)
    % Processes images into a valid input to the program
    % impath indicates the background colour: 'k' is black and 'w' is white
    
    [preProcImage,activityLayer]=ProcImageConstraint(imPath);
    
    [finishedImage]=ResizeImage(preProcImage);
    [activityLayer]=ResizeImage(activityLayer);
    [finishedImage]=AddZeroLayer(finishedImage,cellLevels);
    cellLevels=size(finishedImage);
    [constrVal]=ArrayToConstraint(finishedImage,activityLayer);
end

function [constrVal]=ArrayToConstraint(finishedImage,activityLayer)
    
    [indConst]=find(activityLayer);
    
    constrVal={indConst',[finishedImage(indConst)]'};
    
    
end

function [preProcImage,activityLayer]=ProcImageConstraint(imPath)
    % Load Image and reduce it to an averaged double array from 0 to 1
    imPath=MakePathCompliant(imPath);
    preProcImage=imread(imPath);
    imClass=class(preProcImage);
    numBit=str2num(regexprep(imClass,'uint',''));
    
    preProcImage=double(preProcImage);
    coarseConstraint=preProcImage(:,:,1); % red
    fineConstraint=preProcImage(:,:,2); % green
    activityLayer=preProcImage(:,:,3)==0; % Blue
    
    lvlMax=2^numBit-1;
    lvlMin=0;
    
    if lvlMin~=lvlMax
        preProcImage=(coarseConstraint*lvlMax+fineConstraint)/lvlMax^2;
    else
        preProcImage=zeros(size(preProcImage));
    end
    
end

function [finishedImage]=ResizeImage(procImage)
    % Resizes the image to the right size for the execution of the program
    
    %finishedImage=TrimZeros(procImage);
    %finishedImage=SquareArray(procImage);
    
    finishedImage=procImage';
    finishedImage=finishedImage(:,end:-1:1);
end

function [array]=AddZeroLayer(array,cellLevels)
    % Adds n layer of zeros to every side of the array
    
    sizArray=size(array);
    nPad=(cellLevels-sizArray);
    nPad(nPad<0)=0;
    nPad=round(nPad/2);
    
    nDim=length(sizArray);
    for ii=1:nDim
        sizArray=size(array);
        sizePads=sizArray;
        sizePads(ii)=nPad(ii);
        pad=zeros(sizePads);
        array=cat(ii,pad,array);
        array=cat(ii,array,pad);
    end
end


%% Design Variable cases

function [population]=DesignVariableConsCaller(constrName,constrVal,paroptim,population,varargin)
    
    switch constrName
        case 'MeanVolFrac'
            [population]=MeanVolumeFraction(constrVal,paroptim,population,varargin{1});
        case 'MinSumVolFrac'
            [population]=MinSumVolumeFraction(constrVal,paroptim,population,varargin{1});
        case 'Naca0012'
            [constrVal]=NacaOuterLimit0012(varargin{1},paroptim);
            [population]=LocalVolumeFraction_min(constrVal,population);
            
        case 'LocalVolFrac_min'
            [population]=LocalVolumeFraction_min(constrVal,population);
        case 'LocalVolFrac_equal'
            [population]=LocalVolumeFraction_equal(constrVal,population);
        case 'LocalVolFrac_max'
            [population]=LocalVolumeFraction_max(constrVal,population);
        case ' '
            
        otherwise
            error('Design Variable Constraint Not Recognised')
    end
    
    
    
end

function [population]=MeanVolumeFraction(constrVal,paroptim,population,baseGrid)
    
    varExtract={'desVarRange'};
    [desVarRange]=ExtractVariables(varExtract,paroptim);
    cellCentred=CellCentreGridInformation(baseGrid);
    isActive=logical([cellCentred(:).isactive]);
    volVec=[cellCentred(isActive).isactive];
    totvol=sum(volVec);
    
    for ii=1:length(population)
        
        fillStart=population(ii).fill;
        
        meanFill=(sum(fillStart.*volVec)/totvol);
        ratio=constrVal/meanFill;
        if ratio<=1
            population(ii).fill=fillStart*ratio;
        else
            maxFill=max(fillStart);
            if maxFill*ratio<=max(desVarRange)
                population(ii).fill=fillStart*ratio;
            else
                
                [population(ii).fill,population(ii).constraint]=...
                    IterativeMeanFill(fillStart,desVarRange,constrVal,volVec,totvol);
                
            end
        end
        
        
    end
    
    
    
end

function [constrVal]=NacaOuterLimit0012(gridrefined,paramoptim)
    
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
    figure
    hold on
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

function [population]=LocalVolumeFraction_min(constrVal,population)
    
    
    indConstr=constrVal{1};
    valConstr=constrVal{2};
    for ii=1:length(population)
        
        fillStart=population(ii).fill;
        
        nonSatConstr=fillStart(indConstr)<valConstr;
        
        fillStart(indConstr(nonSatConstr))=valConstr(nonSatConstr);
        population(ii).fill=fillStart;
    end
    
    
    
end

function [population]=LocalVolumeFraction_max(constrVal,population)
    
    
    indConstr=constrVal{1};
    valConstr=constrVal{2};
    for ii=1:length(population)
        
        fillStart=population(ii).fill;
        
        nonSatConstr=fillStart(indConstr)>valConstr;
        
        fillStart(indConstr(nonSatConstr))=valConstr(nonSatConstr);
        population(ii).fill=fillStart;
    end
end

function [population]=LocalVolumeFraction_equal(constrVal,population)
    
    
    indConstr=constrVal{1};
    valConstr=constrVal{2};
    for ii=1:length(population)
        
        fillStart=population(ii).fill;
        
        nonSatConstr=fillStart(indConstr)~=valConstr;
        
        fillStart(indConstr(nonSatConstr))=valConstr(nonSatConstr);
        population(ii).fill=fillStart;
    end
end

function [fill,isConstr]=IterativeMeanFill(fill,desVarRange,constrVal,volVec,totvol)
    isConstr=true;
    ratio=constrVal/(sum(fillStart.*volVec)/totvol);
    kk=0;
    n=length(fill);
    while ratio~=constrVal && kk<=n+1;
        maxFill=max(desVarRange);
        
        fillBound=((fill*ratio)>=maxFill);
        fill(fillBound)=maxFill;
        fill(~fillBound)=fill(~fillBound)*ratio;
        ratio=constrVal/(sum(fill.*volVec)/totvol);
        kk=kk+1;
    end
    if kk>n+1
        isConstr=false;
    end
end

function [population]=MinSumVolumeFraction(constrVal,paroptim,population,baseGrid)
    
    varExtract={'desVarRange'};
    [desVarRange]=ExtractVariables(varExtract,paroptim);
    cellCentred=CellCentreGridInformation(baseGrid);
    isActive=logical([cellCentred(:).isactive]);
    volVec=[cellCentred(isActive).isactive];
    totvol=sum(volVec);
    
    for ii=1:length(population)
        
        fillStart=population(ii).fill;
        
        sumFill=sum(fillStart.*volVec)/totvol;
        ratio=constrVal/sumFill;
        if ~isfinite(ratio)
            population(ii).fill=fillStart;
            population(ii).constraint=false;
        elseif ratio<=1
            population(ii).fill=fillStart;
        elseif ratio>1
            maxFill=max(fillStart);
            if maxFill*ratio<=max(desVarRange)
                population(ii).fill=fillStart*ratio;
            else
                
                [population(ii).fill,population(ii).constraint]=...
                    IterativeMinFill(fillStart,desVarRange,constrVal,volVec,totvol);
                
            end
            
        end
        
        
    end
    
    
    
end

function [fill,isConstr]=IterativeMinFill(fill,desVarRange,constrVal,volVec,totvol)
    isConstr=true;
    ratio=constrVal/sum(fill);
    kk=0;
    n=length(fill);
    while ratio>1 && kk<=n+1;
        maxFill=max(desVarRange);
        
        fillBound=((fill*ratio)>=maxFill);
        fill(fillBound)=maxFill;
        fill(~fillBound)=fill(~fillBound)*ratio;
        ratio=constrVal/sum(fill.*volVec)*totvol;
        kk=kk+1;
    end
    if kk>n+1
        isConstr=false;
    end
end

%% Results cases

function [population]=ResultVariableConsCaller(constrName,constrVal,paroptim,population,varargin)
    
    switch constrName
        case 'AeroResidual'
            population=CheckAerodynamicResidual(constrVal,population);
        case 'AeroResidualBarrier'
            population=BarrierAerodynamicResidual(constrVal,population);
        case 'AeroLift'
            
        case 'Volume'
            
        case 'SnakResBarrier'
            constrVal=ExtractVariables({'convLevel'},paroptim.parametrisation);
            constrVal=[constrVal*100,1];
            population=BarrierSnaxelVolResidual(constrVal,population);
        case ' '
            
        otherwise
            error('Design Variable Constraint Not Recognised')
    end
    
    
    
end

function population=CheckAerodynamicResidual(constrVal,population)
    
    for ii=1:length(population)
        
        if population(ii).constraint
            constrViolation= (population(ii).additional.res>constrVal) && ...
                (population(ii).additional.res~=0);
            population(ii).constraint=~constrViolation;
        end
        
    end
    
    
end

function population=BarrierAerodynamicResidual(constrVal,population)
    
    for ii=1:length(population)
        
        if population(ii).constraint>0
            if population(ii).additional.res<constrVal(1) ... % Constraint fully satisfied
                    || population(ii).additional.res==0
                population(ii).constraint=1;
                
            elseif population(ii).additional.res>constrVal(2) % Constraint fully violated
                population(ii).constraint=0;
                
            else
                population(ii).constraint=1-(0.5+0.5*tanh(3/(-(constrVal(1)-constrVal(2))/2)...
                    *population(ii).additional.res));
            end
        end
        
    end
    
    
end

function population=BarrierSnaxelVolResidual(constrVal,population)
    
    for ii=1:length(population)
        
        if population(ii).constraint>0
            if population(ii).additional.snaxelVolRes<constrVal(1)
                population(ii).constraint=1;
                
            elseif population(ii).additional.snaxelVolRes>constrVal(2) % Constraint fully violated
                population(ii).constraint=0;
                
            else
                population(ii).constraint=1-(0.5+0.5*tanh(3/(-(log10(constrVal(1))-log10(constrVal(2)))/2)...
                    *log10(population(ii).additional.snaxelVolRes)));
            end
        end
        
    end
    
    
end


