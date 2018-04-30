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
        'nSurfPoints','geomStepMag','objFuncCalls','CD0','errorVecMode','location','nTopo'};
    standardIn=cell(size(expectFields));
    [standardIn{:}]=deal(0);
    structBuild=[expectFields;standardIn];
    ASOstruct=struct(structBuild{:});
    
    [subdivLevel1, errorMagnitude, nDV] = ASO.Postproc.subdivData(pathToASO);
    [majorIt, objective, subdivLevel, geomStep] = ASO.Postproc.iterationData(pathToASO);
    [CD0, nFunCall, nSurfPts,errorMode] = ASO.Postproc.ASOData(pathToASO);
    profiles = ASO.Postproc.profileData(pathToASO);
    
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
    ASOstruct.nTopo=numel(profiles);
    
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
    expectFields={'DEIter','majorIt','obj','refLvl','geomErrMag','ASOdesVec',...
        'nSurfPoints','geomStepMag','objFuncCalls','CD0'};
    %              scalar    vec       vec    vec     vec(numLvl)  vec(numLvl)    scalar     vec            scalar     scalar
    
    c=get(axDeOpt,'colororder');
    fColor=[ 0 0 0];
    errVecModes={ASOstruct.errorVecMode};
    lErrVec=unique(errVecModes);
    ASOstructAll=ASOstruct;
    for iii=1:numel(lErrVec)
        ASOstruct=ASOstructAll(~cellfun(@isempty,regexp(errVecModes,lErrVec{iii})));
        
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
        [summaryStruct(iii)]=PlotASOPerformance_DataExtraction(ASOstruct,lErrVec{iii});
    end
    
    % Generate figures
    if ~isaxDef
        figNames={'ASO Performance','ASO Performance Normalised','ASO desvar to surface'};
        tfunc=@(x)numel(unique(x));
        runExtraFig=any(cellfun(tfunc,{summaryStruct.subDivLevel}));
        if runExtraFig
            
            
            xStr={'subDivLevel','nDesVarperBody'};
            for jj=1:numel(xStr)
                kk=numel(figNames);
                figNames(kk+(1:2))={['ASO Performance (',xStr{jj},')'],...
                    ['ASO Performance Normalised (',xStr{jj},')']};
            end
        end
        
        kk=0;
        for ii=1:numel(figNames)
            h(ii)=figure('Name',figNames{ii});
            
            for jj=1:4
                kk=kk+1;
                ax(kk)=subplot(2,2,jj);
            end
        end
    end
    
    % Plot Data
    fColor=[ 0 0 0];
    fMarker='.';
    for iii=1:numel(summaryStruct)
        if numel(lErrVec)>1
            fColor=c(iii,:);
        end
        PlotASOPerformance_fig1(ax(1:4),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig2(ax(5:8),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        PlotASOPerformance_fig3(ax(9:12),summaryStruct(iii),lErrVec{iii},fMarker,fColor)
        
        % Variable subdiv plots
        kk=12;
        for jj=1:numel(xStr)
            x=summaryStruct(iii).(xStr{jj});
            if runExtraFig
                PlotASOPerformance_asoperffigx(ax(kk+(1:4)),x,xStr{jj},...
                    summaryStruct(iii),lErrVec{iii},fMarker,fColor)
                PlotASOPerformance_asoperffignormx(ax(kk+(5:8)),x,xStr{jj},...
                    summaryStruct(iii),lErrVec{iii},fMarker,fColor)
                kk=kk+8;
            end
        end
        
    end
    
    for ii=[2 3 5 6 7 11 12]
        ax(ii).XScale='log';
    end
    % err mag plots
    AverageAllLines(ax(2),[1e-3 1e-2 1e-1 1 10])
    AverageAllLines(ax(3),[1e-3 1e-2 1e-1 1 10])
    AverageAllLines(ax(6),[1e-3 1e-2 1e-1 1 10])
    AverageAllLines(ax(7),[1e-3 1e-2 1e-1 1 10])
    
    
    AverageAllLines(ax(10),{750:250:2000 , 50:50:250})
    % sub div
    AverageAllLines(ax(14),[1.5:4.5])
    AverageAllLines(ax(15),[1.5:4.5])
    AverageAllLines(ax(18),[1.5:4.5])
    AverageAllLines(ax(19),[1.5:4.5])
    AverageAllLines(ax(23),[50:50:200])
    AverageAllLines(ax(27),[50:50:200])
    AverageAllLines(ax(24),[50:50:200])
    AddNumbersInRange(ax(3),0)
    AddNumbersInRange(ax(11),0)
    AddNumbersInRange(ax(12),0)
    legend(ax(1),findobj(ax(1),'type','line'))
    PlotSummaryStruct(summaryStruct)
    
end

function [summaryStruct]=PlotASOPerformance_DataExtraction(ASOstruct,lErrVec)
    
    if true
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
    nDesVarperBody=[];
    changeObjInitial=[];
    changeObjInitialRat=[];
    changeObjCD0=[];
    changeObjFinalCD0=[];
    geomStepMag=[];
    nSurfPoints=[];
    
    if true
        d2=who;
        keepVar=true(size(d2));
        for ii=1:numel(d)
            keepVar=keepVar & cellfun(@(x)~strcmp(x,d{ii}),d2);
        end
        keepCell=[d2(keepVar);'errorVecMode']';
        standardIn=cell(size(keepCell));
        
        structBuild=[keepCell;standardIn];
        summaryStruct=struct(structBuild{:});
    end
    
    deltaPercent=0.9;
    
    for ii=1:numel(ASOstruct)
        jjStart=numel(errMagVec);
        errMagVec=[errMagVec,ASOstruct(ii).geomErrMag];
        nDesVar=[nDesVar,ASOstruct(ii).ASOdesVec];
        nDesVarperBody=[nDesVarperBody,ASOstruct(ii).ASOdesVec/ASOstruct(ii).nTopo];
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
        summaryStruct.(fieldsSum{ii})=eval(fieldsSum{ii});
    end
    summaryStruct.(fieldsSum{end})=lErrVec;
    
end

function []=PlotASOPerformance_fig1(ax,summarrystruct,caseName,fMarker,fColor)
    
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),summarrystruct.subDivLevel,summarrystruct.errMagVec,fMarker,'color',...
        fColor,'DisplayName',caseName)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),summarrystruct.errMagVec,summarrystruct.changeObj,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),summarrystruct.errMagVec,summarrystruct.changeObjInitial,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),summarrystruct.nDesVar,summarrystruct.nIter95Delta,fMarker,'color',fColor)
    
    
    ax(1).XLabel.String='subdivision level';
    ax(1).YLabel.String='subdiv error mag';
    ax(2).XLabel.String='subdiv error mag';
    ax(2).YLabel.String='Obj Value change';
    ax(3).XLabel.String='subdiv error mag';
    ax(3).YLabel.String='Obj Value change over majorIt 1';
    ax(4).XLabel.String='# Design Variables';
    ax(4).YLabel.String='nIter to 90% obj change';
    
end

function []=PlotASOPerformance_fig2(ax,summarrystruct,caseName,fMarker,fColor)
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),summarrystruct.errMagVec,summarrystruct.nDesVar,fMarker,'color',fColor)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),summarrystruct.errMagVec,summarrystruct.changeObjRat,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),summarrystruct.errMagVec,summarrystruct.changeObjInitialRat,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),summarrystruct.nDesVar,summarrystruct.nIter95DeltaRat,fMarker,'color',fColor)
    
    
    
    ax(1).XLabel.String='subdiv error mag';
    ax(1).YLabel.String='# Design Variables';
    ax(2).XLabel.String='subdiv error mag';
    ax(2).YLabel.String='Obj Value change (% of orig objective)';
    ax(3).XLabel.String='subdiv error mag';
    ax(3).YLabel.String='Obj Value change over majorIt 1 (% total change)';
    ax(4).XLabel.String='# Design Variables';
    ax(4).YLabel.String='nIter to 90% obj change (% total #iter)';
    
end

function []=PlotASOPerformance_fig3(ax,summarrystruct,caseName,fMarker,fColor)
    
    plot(ax(1),summarrystruct.nSurfPoints,summarrystruct.nDesVar,fMarker,'color',fColor)
    % Magnitude of the geometric step vs number of surf points and # of design variables
    plot3(ax(2),summarrystruct.nSurfPoints,summarrystruct.nDesVar,...
        summarrystruct.geomStepMag,fMarker,'color',fColor)
    plot(ax(3),summarrystruct.errMagVec,summarrystruct.changeObjCD0,fMarker,'color',fColor)
    plot(ax(4),summarrystruct.errMagVec,summarrystruct.changeObjFinalCD0,fMarker,'color',fColor)
    % Add lines and numbers indicating overUnder
    
    
    ax(1).XLabel.String='# surf points';
    ax(1).YLabel.String='# design variable';
    ax(2).XLabel.String='# surf points';
    ax(2).YLabel.String='# design variable';
    ax(2).ZLabel.String='Size of Geometric Step';
    ax(3).XLabel.String='subdiv error mag';
    ax(3).YLabel.String='CD0-obj(1)';
    ax(4).XLabel.String='subdiv error mag';
    ax(4).YLabel.String='CD0-obj(end)';
end


function []=PlotASOPerformance_asoperffigx(ax,x,xStr,summarrystruct,caseName,fMarker,fColor)
    
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),x,summarrystruct.errMagVec,fMarker,'color',...
        fColor,'DisplayName',caseName)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),x,summarrystruct.changeObj,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),x,summarrystruct.changeObjInitial,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),x,summarrystruct.nIter95Delta,fMarker,'color',fColor)
    
    for ii=1:4
        [ax(ii).XLabel.String]=(xStr);
    end
    ax(1).YLabel.String='subdiv error mag';
    ax(2).YLabel.String='Obj Value change';
    ax(3).YLabel.String='Obj Value change over majorIt 1';
    ax(4).YLabel.String='nIter to 90% obj change';
    
end

function []=PlotASOPerformance_asoperffignormx(ax,x,xStr,summarrystruct,caseName,fMarker,fColor)
    
    % Error Magnitude (and later subdivision level) vs
    % 1) geometric error Magnitude (cloud)
    plot(ax(1),x,summarrystruct.nDesVar,fMarker,'color',fColor)
    % 2) Change in objective thanks to ASO (cloud)
    plot(ax(2),x,summarrystruct.changeObjRat,fMarker,'color',fColor)
    % 3) Initial Objective Delta (should be 0 except for error = 'none')
    plot(ax(3),x,summarrystruct.changeObjInitialRat,fMarker,'color',fColor)
    % 4) Change of objective other the first iteration (normalised by total Delta)
    plot(ax(4),x,summarrystruct.nIter95DeltaRat,fMarker,'color',fColor)
    
    for ii=1:4
        [ax(ii).XLabel.String]=(xStr);
    end
    ax(1).YLabel.String='# Design Variables';
    ax(2).YLabel.String='Obj Value change (% of orig objective)';
    ax(3).YLabel.String='Obj Value change over majorIt 1 (% total change)';
    ax(4).YLabel.String='nIter to 90% obj change (% total #iter)';
    
end


function []=PlotSummaryStruct(summaryStruct)
    
end


function []=AverageAllLines(ax,ranges)
    hold(ax,'on');
    
    l=findobj(ax,'type','line');
    for ii=1:numel(l)
        if isempty(l(ii).ZData)
            [laverage]=AverageLines(ax,l(ii).XData,l(ii).YData,ranges);
        else
            [laverage]=AverageLines3D(ax,l(ii).XData,l(ii).YData,...
                l(ii).ZData,ranges{1},ranges{2});
        end
        [laverage.LineStyle]=deal('--');
        laverage(1).Marker='*';
        laverage(2).Marker='d';
        [laverage.Color]=deal(l(ii).Color);
    end
end

function [l]=AverageLines(ax,x,y,ranges)
    hold(ax,'on');
    cellRanges=cell(numel(ranges)+1,1);
    ranges=[-inf,ranges,inf];
    xPlotMean=[];
    yPlotMean=[];
    xPlotMedian=[];
    yPlotMedian=[];
    for ii=1:numel(cellRanges)
        cellRanges{ii}=find(x<ranges(ii+1) & x>=ranges(ii));
        if numel(cellRanges{ii})>0
            xPlotMean=[xPlotMean,mean(x(cellRanges{ii}))];
            yPlotMean=[yPlotMean,mean(y(cellRanges{ii}))];
            xPlotMedian=[xPlotMedian,median(x(cellRanges{ii}))];
            yPlotMedian=[yPlotMedian,median(y(cellRanges{ii}))];
        end
    end
    
    l=plot(ax,xPlotMean,yPlotMean,xPlotMedian,yPlotMedian);
end

function [l]=AverageLines3D(ax,x,y,z,rangesX,rangesY)
    hold(ax,'on');
    cellRanges=cell(numel(rangesX)+1,numel(rangesY)+1);
    rangesX=[-inf,rangesX,inf];
    rangesY=[-inf,rangesY,inf];
    
    xPlotMean=zeros(numel(rangesX)-1,numel(rangesY)-1);
    yPlotMean=zeros(numel(rangesX)-1,numel(rangesY)-1);
    xPlotMedian=zeros(numel(rangesX)-1,numel(rangesY)-1);
    yPlotMedian=zeros(numel(rangesX)-1,numel(rangesY)-1);
    zPlotMean=zeros(numel(rangesX)-1,numel(rangesY)-1);
    zPlotMedian=zeros(numel(rangesX)-1,numel(rangesY)-1);
    
    for ii=1:(numel(rangesX)-1)
        for jj=1:(numel(rangesY)-1)
            cellRanges{ii,jj}=find(x<rangesX(ii+1) & x>=rangesX(ii) ...
                & y<rangesY(jj+1) & x>=rangesX(jj));
            if numel(cellRanges{ii})>0
                xPlotMean(ii,jj)=mean(x(cellRanges{ii,jj}));
                yPlotMean(ii,jj)=mean(y(cellRanges{ii,jj}));
                zPlotMean(ii,jj)=mean(z(cellRanges{ii,jj}));
                
                xPlotMedian(ii,jj)=median(x(cellRanges{ii,jj}));
                yPlotMedian(ii,jj)=median(y(cellRanges{ii,jj}));
                zPlotMedian(ii,jj)=median(z(cellRanges{ii,jj}));
            end
        end
    end
    
    %l=plot(ax,xPlotMean,yPlotMean,xPlotMedian,yPlotMedian);
    l(1)=mesh(ax,xPlotMean,yPlotMean,zPlotMean);
    l(2)=mesh(ax,xPlotMedian,yPlotMedian,zPlotMedian);
end

function []=AddNumbersInRange(ax,overUnder)
    hold(ax,'on')
    flag=true;
    l=findobj(ax,'type','line');
    
    y=[l.YData];
    over=sum(y>overUnder);
    under=numel(y)-over;
    
    boxAx=axis(ax);
    boxAx(1:2)=[min([l.XData]),max([l.XData])];
    if strcmp(ax.XScale,'linear')
        boxAx(1:2)=boxAx(1:2)+[-0.1 0.1]*(boxAx(2)-boxAx(1));
    else
        boxAx(1:2)=10.^(log10(boxAx(1:2))+([-0.1 0.1]*log10(boxAx(2)-boxAx(1))));
    end
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
        plot(ax,[min([l.XData]),max([l.XData])],[1 1]*overUnder,'--','color',[0.4 0.4 0.4])
    end
    axis(ax,boxAx)
end



