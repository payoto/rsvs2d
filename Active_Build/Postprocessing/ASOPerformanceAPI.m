function [ASOstruct,h]=ASOPerformanceAPI(optIn, ASOiters,dirSave)
    % prepare for different cases
    runASOextract=1;
    if isstruct(optIn)
        try
            optimstruct=optIn.optimstruct;
            paramoptim=optIn.paramoptim;
            
        catch
            runASOextract=0;
            ASOstruct=optIn;
        end
    elseif ischar(optIn)
        optPath=FindDir(optIn,'OptimRes',0);
        paramPath=FindDir(optIn,'FinalParam',0);
        optimstruct=load(optPath{1});
        optimstruct=optimstruct.optimstruct;
        paramoptim=load(paramPath{1});
        paramoptim=paramoptim.paramoptim;
        if nargin<3
           dirSave=optIn; 
        end
    elseif iscell(optIn)
        
        paramPath=FindDir(optIn{1},'FinalParam',0);
        
        paramoptim=load(paramPath{1});
        paramoptim=paramoptim.paramoptim;
        
        for ii=1:numel(optIn)
            optPath=FindDir(optIn{ii},'OptimRes',0);
            loadDat=load(optPath{1});
            if exist('optimstruct','var')
                optimstruct=[optimstruct,loadDat.optimstruct];
            else
                optimstruct=loadDat.optimstruct;
            end
        end
        
    else
        error('unknown input argument type')
    end
    % Do stuff
    knownOptim=0;
    
    kk=1;
    ll=0;
    if runASOextract
        
        varExtract={'direction','knownOptim','objectiveName'};
        [direction,knownOptim,objectiveName]=ExtractVariables(varExtract,paramoptim);
        % ---------
        % This will cause issues when using many optimstructs stuck
        % together
        rmPop=false(size(optimstruct));
        for ii=1:numel(optimstruct)
            rmPop(ii)=isempty(optimstruct(ii).population(1).objective);
        end
        optimstruct=optimstruct(~rmPop);
        % ---------
        for ii=ASOiters
            for jj=1:numel(optimstruct(ii).population)
                try
                    ASOstruct(kk)=ASOInterface(optimstruct(ii).population(jj).location);
                    ASOstruct(kk).DEIter=ii-1;
                    kk=kk+1;
                catch MEid
                    MEid.getReport
                    disp([int2str(ii),':',int2str(jj),' failed'])
                    ll=ll+1;
                end
                
            end
        end
        disp([int2str(ll),' failures to ',int2str(kk),' success'])
        h=OptimHistory(0,optimstruct,knownOptim,1000,'min');
        ax=findobj(h(1),'type','axes');
        if exist('dirSave','var')
            
            save([dirSave,filesep,'ASOperformance.mat'],'ASOstruct');
            assignin('base','ASOstruct',ASOstruct);
        end
    else
        h=figure;
        ax(2)=axes;
        hold on
    end
    [h2]=PlotASOPerformance(ASOstruct,ax(2));
    h=[h,h2];
    
    % Save Data
    if exist('dirSave','var')
        
        save([dirSave,filesep,'ASOperformance.mat'],'ASOstruct');
        hgsave(h,[dirSave,filesep,'ASOperformanceFigures.mat'])
        
    end
    
end

function [ASOstruct]=ASOInterface(pathToASO)
    expectFields={'DEIter','majorIt','obj','refLvl','geomErrMag','ASOdesVec',...
        'nSurfPoints','geomStepMag','objFuncCalls','CD0'};
    standardIn=cell(size(expectFields));
    [standardIn{:}]=deal(0);
    structBuild=[expectFields;standardIn];
    ASOstruct=struct(structBuild{:});
    
    [subdivLevel1, errorMagnitude, nDV] = ASO.Postproc.subdivData(pathToASO);
    [majorIt, objective, subdivLevel, geomStep] = ASO.Postproc.iterationData(pathToASO);
    [CD0, nFunCall, nSurfPts] = ASO.Postproc.ASOData(pathToASO);
    
    ASOstruct.DEIter=1;
    ASOstruct.majorIt=majorIt;
    ASOstruct.obj=objective;
    ASOstruct.refLvl=unique(subdivLevel);
    ASOstruct.geomErrMag=errorMagnitude(FindObjNum([],unique(subdivLevel),(subdivLevel1)));
    ASOstruct.ASOdesVec=nDV(FindObjNum([],unique(subdivLevel),(subdivLevel1)))';
    ASOstruct.nSurfPoints=nSurfPts;
    ASOstruct.geomStepMag=geomStep;
    ASOstruct.objFuncCalls=nFunCall;
    ASOstruct.CD0=CD0;
    if isempty(majorIt)
       error('No ASO was run') 
    end
    
    if ~isfinite(objective(1))
       error('ASO failed') 
    end
    
end

function [h]=PlotASOPerformance(ASOstruct,axDeOpt)
    
    
    fieldsASO=fieldnames(ASOstruct);
    expectFields={'DEIter','majorIt','obj','refLvl','geomErrMag','ASOdesVec','nSurfPoints','geomStepMag','objFuncCalls','CD0'};
    %              scalar    vec       vec    vec     vec(numLvl)  vec(numLvl)    scalar     vec            scalar     scalar
    
    % Iteration plots overlaid on DEresults
    for ii=1:numel(ASOstruct)
        plot(axDeOpt,[ASOstruct(ii).DEIter,ASOstruct(ii).DEIter+1],...
            ASOstruct(ii).obj([1,end]),'k--');
        % Possibly change this to have a changing color with better
        % improvements
    end
    
    % Data rearrangement
    errMagVec=[];
    subDivLevel=[];
    changeObj=[];
    changeObjRat=[];
    nIter95Delta=[];
    nIter95DeltaRat=[];
    nDesVar=[];
    changeObjInitial=[];
    changeObjInitialRat=[];
    changeObjCD0=[];
    geomStepMag=[];
    deltaPercent=0.9;
    for ii=1:numel(ASOstruct)
        jjStart=numel(errMagVec);
        errMagVec=[errMagVec,ASOstruct(ii).geomErrMag];
        nDesVar=[nDesVar,ASOstruct(ii).ASOdesVec];
        subDivLevel=[subDivLevel,unique(ASOstruct(ii).refLvl)];
        
        if numel(errMagVec)~=numel(subDivLevel)
            error('subDivLevel size does not match number of errMagnitudes provided')
        end
        for jj=1:(numel(errMagVec)-jjStart)
            currList=find(ASOstruct(ii).refLvl==subDivLevel(jj+jjStart));
            if numel(currList)==1
                currList(2)=currList(1);
            end
            changeObj(jjStart+jj)=min(ASOstruct(ii).obj(currList))...
                -ASOstruct(ii).obj(currList(1));
            changeObjRat(jjStart+jj)=(min(ASOstruct(ii).obj(currList))...
                -ASOstruct(ii).obj(currList(1)))/ASOstruct(ii).obj(currList(1));
            changeObjInitial(jjStart+jj)=(ASOstruct(ii).obj(currList(2))...
                -ASOstruct(ii).obj(currList(1)));
            changeObjInitialRat(jjStart+jj)=(ASOstruct(ii).obj(currList(2))...
                -ASOstruct(ii).obj(currList(1)))/(ASOstruct(ii).obj(currList(end))...
                -ASOstruct(ii).obj(currList(1)));
            
            nIter95Delta(jjStart+jj)=min(find((ASOstruct(ii).obj(1)+changeObj(jjStart+jj)*deltaPercent)>=ASOstruct(ii).obj(currList)));
            nIter95DeltaRat(jjStart+jj)=min(find((ASOstruct(ii).obj(1)+changeObj(jjStart+jj)*...
                deltaPercent)>=ASOstruct(ii).obj(currList)))/numel(currList);
            nSurfPoints(jj+jjStart)=ASOstruct(ii).nSurfPoints;
            
            
            changeObjCD0(jj+jjStart)=ASOstruct(ii).obj(currList(1))-ASOstruct(ii).CD0;
            geomStepMag(jj+jjStart)=ASOstruct(ii).geomStepMag(currList(end));
        end
        
    end
    
    h=figure('Name','ASO Performance');
    for ii=1:4
        ax(ii)=subplot(2,2,ii);
    end
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),subDivLevel,errMagVec,'ko')
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),errMagVec,changeObj,'ko')
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),errMagVec,changeObjInitial,'ko')
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),nDesVar,nIter95Delta,'ko')
    
    h(2)=figure('Name','ASO Performance Normalised');
    for ii=5:8
        ax(ii)=subplot(2,2,ii);
    end
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(5),subDivLevel,errMagVec,'ko')
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(6),errMagVec,changeObjRat,'ko')
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(7),errMagVec,changeObjInitialRat,'ko')
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(8),nDesVar,nIter95DeltaRat,'ko')
    
    % #iterations needed to achieve 95% of the improvement vs Number of Desvar
    
    h(3)=figure('name','ASO desvar to surface');
    for ii=9:12
        ax(ii)=subplot(2,2,ii);
    end
    plot(ax(9),nSurfPoints,nDesVar,'ko')
    % Magnitude of the geometric step vs number of surf points and # of design variables
    plot3(ax(10),nSurfPoints,nDesVar,geomStepMag,'ko')
    plot(ax(11),errMagVec,changeObjCD0,'ko')
end