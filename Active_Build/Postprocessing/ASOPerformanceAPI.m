function [ASOstruct,h]=ASOPerformanceAPI(optIn, ASOiters,dirSave,nameRun)
    % prepare for different cases
    runASOextract=1;
    if nargin<4
        nameRun='';
    end
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
        if nargin<4
            nameRun=regexprep(optIn,'^.*Dir[0-9,\-,_]*T[0-9]*_','');
        end
    elseif iscell(optIn)
        
        paramPath=FindDir(optIn{1},'FinalParam',0);
        
        paramoptim=load(paramPath{1});
        paramoptim=paramoptim.paramoptim;
        
        for ii=1:numel(optIn)
            optPath=FindDir(optIn{ii},'OptimRes',0);
            loadDat=load(optPath{1});
            if exist('optimstruct','var')
                optimstruct=[optimstruct,loadDat.optimstruct(ASOiters)];
            else
                optimstruct=loadDat.optimstruct;
            end
            
        end
        ASOiters=1:numel(optimstruct);
        
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
                    ASOstruct(kk).location=optimstruct(ii).population(jj).location;
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
            
            save([dirSave,filesep,'ASOperformance',nameRun,'.mat'],'ASOstruct');
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
        
        save([dirSave,filesep,'ASOperformance',nameRun,'.mat'],'ASOstruct');
        for ii=1:numel(h)
            h(ii).Name= [h(ii).Name,' - ',nameRun];
        end
        hgsave(h,[dirSave,filesep,'ASOperformance',nameRun,'.fig'])
        
    end
    
end

function [ASOstruct]=ASOInterface(pathToASO)
    expectFields={'DEIter','majorIt','obj','refLvl','geomErrMag','ASOdesVec',...
        'nSurfPoints','geomStepMag','objFuncCalls','CD0','errorVecMode','location'};
    standardIn=cell(size(expectFields));
    [standardIn{:}]=deal(0);
    structBuild=[expectFields;standardIn];
    ASOstruct=struct(structBuild{:});
    
    [subdivLevel1, errorMagnitude, nDV] = ASO.Postproc.subdivData(pathToASO);
    [majorIt, objective, subdivLevel, geomStep] = ASO.Postproc.iterationData(pathToASO);
    [CD0, nFunCall, nSurfPts,errorMode] = ASO.Postproc.ASOData(pathToASO);
    
    ASOstruct.DEIter=1;
    ASOstruct.majorIt=majorIt;
    ASOstruct.obj=objective;
    ASOstruct.refLvl=(subdivLevel);
    ASOstruct.geomErrMag=errorMagnitude(FindObjNum([],unique(subdivLevel),(subdivLevel1)));
    ASOstruct.ASOdesVec=nDV(FindObjNum([],unique(subdivLevel),(subdivLevel1)))';
    ASOstruct.nSurfPoints=nSurfPts;
    ASOstruct.geomStepMag=geomStep;
    ASOstruct.objFuncCalls=nFunCall;
    ASOstruct.CD0=CD0;
    ASOstruct.errorVecMode=errorMode;
    if isempty(majorIt)
        error('No ASO was run')
    end
    
    if ~isfinite(objective(1))
        error('ASO failed')
    end
    
end

function [h,ax]=PlotASOPerformance(ASOstruct,axDeOpt,axOther)
    isaxDef=false;
    if nargin==3
        isaxDef=true;
        ax=axOther;
    end
    fieldsASO=fieldnames(ASOstruct);
    expectFields={'DEIter','majorIt','obj','refLvl','geomErrMag','ASOdesVec','nSurfPoints','geomStepMag','objFuncCalls','CD0'};
    %              scalar    vec       vec    vec     vec(numLvl)  vec(numLvl)    scalar     vec            scalar     scalar
    
    c=get(axDeOpt,'colororder');
    fColor=[ 0 0 0];
    errVecModes={ASOstruct.errorVecMode};
    lErrVec=unique(errVecModes);
    ASOstructAll=ASOstruct;
    for iii=1:numel(lErrVec)
        ASOstruct=ASOstructAll(~cellfun(@isempty,regexp(errVecModes,lErrVec{iii})));
        if iii>1
            isaxDef=true;
            for ii=1:numel(ax)
                hold(ax(ii),'on')
            end
        end
        if numel(lErrVec)>1
            fColor=c(iii,:);
        end
        % Iteration plots overlaid on DEresults
        for ii=1:numel(ASOstruct)
            plot(axDeOpt,[ASOstruct(ii).DEIter,ASOstruct(ii).DEIter+1],...
                [ASOstruct(ii).CD0,ASOstruct(ii).obj(end)],'--','color',fColor);
            % Possibly change this to have a changing color with better
            % improvements
        end
        
        % Data rearrangement
        
        if iii==1
            d2=[];
            keepVar=[];
            d=[];
            d=who;
        end
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
        changeObjFinalCD0=[];
        geomStepMag=[];
        nSurfPoints=[];
        
        if iii==1
            d2=who;
            keepVar=true(size(d2));
            for ii=1:numel(d)
                keepVar=keepVar & cellfun(@(x)~strcmp(x,d{ii}),d2);
            end
            keepCell=[d2(keepVar);'errorVecMode']';
            standardIn=cell(size(keepCell));
            
            structBuild=[keepCell;standardIn];
            summaryStruct=struct(structBuild{:});
            summaryStruct=repmat(summaryStruct,size(lErrVec));
        end
        
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
                
                
                changeObjCD0(jj+jjStart)=ASOstruct(ii).obj(1)-ASOstruct(ii).CD0;
                changeObjFinalCD0(jj+jjStart)=ASOstruct(ii).obj(end)-ASOstruct(ii).CD0;
                geomStepMag(jj+jjStart)=ASOstruct(ii).geomStepMag(currList(end));
            end
            
        end
        fieldsSum=fieldnames(summaryStruct);
        for ii=1:numel(fieldsSum)-1
            summaryStruct(iii).(fieldsSum{ii})=eval(fieldsSum{ii});
        end
        summaryStruct(iii).(fieldsSum{end})=lErrVec{iii};
        
        if ~isaxDef
            h=figure('Name','ASO Performance');
            for ii=1:4
                ax(ii)=subplot(2,2,ii);
            end
        else
            isaxDef=numel(ax)>4;
        end
        % Error Magnitude (and later subdivision level) vs
        % 1) geometric error Magnitude (cloud)
        plot(ax(1),subDivLevel,errMagVec,'o','color',fColor,'DisplayName',lErrVec{iii})
        % 2) Change in objective thanks to ASO (cloud)
        plot(ax(2),errMagVec,changeObj,'o','color',fColor)
        % 3) Initial Objective Delta (should be 0 except for error = 'none')
        plot(ax(3),errMagVec,changeObjInitial,'o','color',fColor)
        % 4) Change of objective other the first iteration (normalised by total Delta)
        plot(ax(4),nDesVar,nIter95Delta,'o','color',fColor)
        if ~isaxDef
            h(2)=figure('Name','ASO Performance Normalised');
            for ii=5:8
                ax(ii)=subplot(2,2,ii-4);
            end
        else
            isaxDef=numel(ax)>8;
        end
        % Error Magnitude (and later subdivision level) vs
        % 1) geometric error Magnitude (cloud)
        plot(ax(5),errMagVec,nDesVar,'o','color',fColor)
        % 2) Change in objective thanks to ASO (cloud)
        plot(ax(6),errMagVec,changeObjRat,'o','color',fColor)
        % 3) Initial Objective Delta (should be 0 except for error = 'none')
        plot(ax(7),errMagVec,changeObjInitialRat,'o','color',fColor)
        % 4) Change of objective other the first iteration (normalised by total Delta)
        plot(ax(8),nDesVar,nIter95DeltaRat,'o','color',fColor)
        
        
        
        ax(1).XLabel.String='subdivision level';
        ax(1).YLabel.String='subdiv error mag';
        ax(2).XLabel.String='subdiv error mag';
        ax(2).YLabel.String='Obj Value change';
        ax(3).XLabel.String='subdiv error mag';
        ax(3).YLabel.String='Obj Value change over majorIt 1';
        ax(4).XLabel.String='# Design Variables';
        ax(4).YLabel.String='nIter to 90% obj change';
        ax(5).XLabel.String='subdiv error mag';
        ax(5).YLabel.String='# Design Variables';
        ax(6).XLabel.String='subdiv error mag';
        ax(6).YLabel.String='Obj Value change (% of orig objective)';
        ax(7).XLabel.String='subdiv error mag';
        ax(7).YLabel.String='Obj Value change over majorIt 1 (% total change)';
        ax(8).XLabel.String='# Design Variables';
        ax(8).YLabel.String='nIter to 90% obj change (% total #iter)';
        % #iterations needed to achieve 95% of the improvement vs Number of Desvar
        if ~isaxDef
            h(3)=figure('name','ASO desvar to surface');
            for ii=9:12
                ax(ii)=subplot(2,2,ii-8);
            end
        else
            isaxDef=numel(ax)>12;
        end
        plot(ax(9),nSurfPoints,nDesVar,'o','color',fColor)
        % Magnitude of the geometric step vs number of surf points and # of design variables
        plot3(ax(10),nSurfPoints,nDesVar,geomStepMag,'o','color',fColor)
        plot(ax(11),errMagVec,changeObjCD0,'o','color',fColor)
        plot(ax(12),errMagVec,changeObjFinalCD0,'o','color',fColor)
        % Add lines and numbers indicating overUnder
        
        
        ax(9).XLabel.String='# surf points';
        ax(9).YLabel.String='# design variable';
        ax(10).XLabel.String='# surf points';
        ax(10).YLabel.String='# design variable';
        ax(10).ZLabel.String='Size of Geometric Step';
        ax(11).XLabel.String='subdiv error mag';
        ax(11).YLabel.String='CD0-obj(1)';
        ax(12).XLabel.String='subdiv error mag';
        ax(12).YLabel.String='CD0-obj(end)';
        
        if numel(unique(subDivLevel))>1
            kk=12;
            x=subDivLevel;
            xStr='Subdivision Level';
            if ~isaxDef
                h(end+1)=figure('Name',['ASO Performance (',xStr,')']);
                
                for ii=1:4
                    ax(kk+ii)=subplot(2,2,ii);
                end
            else
                isaxDef=numel(ax)>(kk+4);
            end
            % Error Magnitude (and later subdivision level) vs
            % 1) geometric error Magnitude (cloud)
            plot(ax(kk+1),x,nDesVar,'o','color',fColor)
            % 2) Change in objective thanks to ASO (cloud)
            plot(ax(kk+2),x,changeObj,'o','color',fColor)
            % 3) Initial Objective Delta (should be 0 except for error = 'none')
            plot(ax(kk+3),x,changeObjInitial,'o','color',fColor)
            % 4) Change of objective other the first iteration (normalised by total Delta)
            plot(ax(kk+4),nDesVar,nIter95Delta,'o','color',fColor)
            if ~isaxDef
                h(end+1)=figure('Name',['ASO Performance Normalised (',xStr,')']);
                for ii=5:8
                    ax(kk+ii)=subplot(2,2,ii-4);
                end
            else
                isaxDef=numel(ax)>(kk+4);
            end
            % Error Magnitude (and later subdivision level) vs
            % 1) geometric error Magnitude (cloud)
            plot(ax(kk+5),x,errMagVec,'o','color',fColor)
            % 2) Change in objective thanks to ASO (cloud)
            plot(ax(kk+6),x,changeObjRat,'o','color',fColor)
            % 3) Initial Objective Delta (should be 0 except for error = 'none')
            plot(ax(kk+7),x,changeObjInitialRat,'o','color',fColor)
            % 4) Change of objective other the first iteration (normalised by total Delta)
            plot(ax(kk+8),nDesVar,nIter95DeltaRat,'o','color',fColor)
            
            ax(kk+1).XLabel.String=xStr;
            ax(kk+1).YLabel.String='design Variables';
            ax(kk+2).XLabel.String=xStr;
            ax(kk+2).YLabel.String='Obj Value change';
            ax(kk+3).XLabel.String=xStr;
            ax(kk+3).YLabel.String='Obj Value change over majorIt 1';
            ax(kk+4).XLabel.String='# Design Variables';
            ax(kk+4).YLabel.String='nIter to 90% obj change';
            ax(kk+5).XLabel.String=xStr;
            ax(kk+5).YLabel.String='subdiv error mag';
            ax(kk+6).XLabel.String=xStr;
            ax(kk+6).YLabel.String='Obj Value change (% of orig objective)';
            ax(kk+7).XLabel.String=xStr;
            ax(kk+7).YLabel.String='Obj Value change over majorIt 1 (% total change)';
            ax(kk+8).XLabel.String='# Design Variables';
            ax(kk+8).YLabel.String='nIter to 90% obj change (% total #iter)';
        end
        
    end
    AddNumbersInRange(ax(11),0)
    AddNumbersInRange(ax(12),0)
    legend(ax(1),findobj(ax(1),'type','line'))
    PlotSummaryStruct(summaryStruct)
end

function []=PlotSummaryStruct(summaryStruct)
    
end

function []=AddNumbersInRange(ax,overUnder)
    hold(ax,'on')
    flag=true;
    l=findobj(ax,'type','line');
    
    y=[l.YData];
    over=sum(y>overUnder);
    under=numel(y)-over;
    
    boxAx=axis(ax);
    
    xDelta=boxAx(2)-0.02*(boxAx(2)-boxAx(1));
    yOver=(boxAx(4)-overUnder)*0.5+overUnder;
    yUnder=(overUnder-boxAx(3))*0.5+boxAx(3);
    if(yUnder)<boxAx(4) && under>0
        text(ax,xDelta,yUnder,int2str(under),'HorizontalAlignment','right','color',[0.4 0.4 0.4]);
    else
        flag=false;
    end
    if(yUnder)>boxAx(3) && over>0
        text(ax,xDelta,yOver,int2str(over),'HorizontalAlignment','right','color',[0.4 0.4 0.4]);
    else
        flag=false;
    end
    if flag
        plot(ax,boxAx(1:2),[1 1]*overUnder,'--','color',[0.4 0.4 0.4])
    end
    axis(ax,boxAx)
end
