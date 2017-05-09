function [resstruct]=PostTreatGroupRun(entryPoint,varargin)
    % This function post treats an entire group of runs. Especially
    % refienement runs
    plotStart=1;
    switch entryPoint
        case 'data'
            resstruct=CollectGroupRunInformation(varargin{:});
            compareopts=varargin{2};
        case 'plot'
            resstruct=varargin{1};
            compareopts=varargin{2};
            if numel(varargin)>2
                plotStart=varargin{3};
            end
    end
    
%     try 
        PlotGroupedInformation(resstruct,compareopts,plotStart)
%     catch ME
%         disp(ME.getReport)
%     end
end

function [resstruct]=CollectGroupRunInformation(rootDir,compareopts)
    
    [childFolder]=ExploreFolderTree({rootDir});
    [optimPaths]=FindOptimResFiles(childFolder);
    %[paramPaths]=FindFinalParamFiles(childFolder);
    
    [resstruct]=TreatPathsIntoStruct(optimPaths,compareopts);
    save([rootDir,filesep,'AllRefineResults.mat'],'resstruct')
end

function []=PlotGroupedInformation(resstruct,compareopts,plotStart)
    
    
    
    if plotStart<2
        for ii=1:numel(resstruct)
            figobj=PlotIndividualRefinement(resstruct(ii),'full');
            hgsave(figobj.h,[resstruct(ii).path,filesep,'FigRefOpt_',...
                resstruct(ii).name{1},'.fig']);
            close(figobj.h)
        end
    end
    optimOpts=[resstruct(:).paramoptim];
    snakeOpts=[resstruct(:).paramsnake];
    [optimOpts,optimList]=IdentifyUniqueOptions(optimOpts);
    [snakeOpts,snakeList]=IdentifyUniqueOptions(snakeOpts);
    [allOpts,allOptsList]=ConcatenateStructures(optimOpts,optimList,...
        snakeOpts,snakeList);
    resDir=resstruct(1).path;
    n=regexp(resDir,filesep);
    resDir=[resDir(1:n(end-1)),'results'];
    mkdir(resDir)
    mkdir([resDir,filesep,'static'])
    
    plotType='ref';
    if plotStart<3
        ExploreParamEffect_SingleStatic(allOptsList,allOpts,resstruct,plotType,...
            [resDir,filesep,'static'])
    end
    mkdir([resDir,filesep,'moving'])
    if plotStart<4
        ExploreParamEffect_SingleMoving(allOptsList,allOpts,resstruct,plotType,...
            [resDir,filesep,'moving'])
    end
    mkdir([resDir,filesep,'carpet'])
    if plotStart<5
        for ii=1:numel(compareopts.paramPair)
            ExtractCarpetParamPair(compareopts.paramPair{ii},allOptsList,allOpts,...
                resstruct,[resDir,filesep,'carpet'])
        end
    end
end

%% Plot Data

function [figobj]=PlotIndividualRefinement(resstruct,typePlot,figobj)
    
    isLoadColor=false;
    if ~exist('figobj','var')
        figobj.h=figure('Name',resstruct.name{1});
        figobj.axh=axes;
        isLoad=false;
    else
        figure(figobj.h)
        axes(figobj.axh)
        isLoadColor=isfield(figobj,'color');
        isLoad=true;
    end
    
    
    hold on
    switch typePlot
        case 'full'
            for ii=1:numel(resstruct.refine)
                
                figobj.l(ii,1)=plot(resstruct.refine(ii).iterstart:resstruct.refine(ii).iterend,...
                    resstruct.refine(ii).obj);
                
                if isLoadColor
                    figobj.l(ii,1).Color=figobj.color;
                end
                figobj.l(ii,1).DisplayName=['Refine Stage ',int2str(resstruct.refine(ii).stage)];
                
                figobj.l(ii,2)=plot(resstruct.refine(ii).iterend,...
                    resstruct.refine(ii).obj(end),'*','Color',figobj.l(ii,1).Color);
                figobj.l(ii,3)=plot([1,resstruct.refine(end).iterend+5],...
                    resstruct.refine(ii).obj(end)*[1 1],'--','Color',figobj.l(ii,1).Color);
                
                figobj.t(ii,1)=text(resstruct.refine(end).iterend+8,...
                    resstruct.refine(ii).obj(end),int2str(resstruct.refine(ii).ndesvar));
                figobj.t(ii,1).VerticalAlignment='middle';
                %figobj.t(ii,1).BackgroundColor=[1 1 1];
                figobj.t(ii,1).Color=figobj.l(ii,1).Color;
                
            end
            
        case 'ref'
            ii=1;
            
            figobj.l(ii,1)=plot([resstruct.refine(:).stage],[resstruct.refine(:).objend]);
            
            if isLoadColor
                figobj.l(ii,1).Color=figobj.color;
            end
            figobj.l(ii,1).DisplayName=['Refine Stage ',int2str(resstruct.refine(ii).stage)];
            figobj.t(ii,1).VerticalAlignment='middle';
            %figobj.t(ii,1).BackgroundColor=[1 1 1];
            figobj.t(ii,1).Color=figobj.l(ii,1).Color;
            
            
    end
    figobj.axh.YScale='log';
    if ~isLoad
        box=axis;
        box(2)=box(2)+10;
        %box(3)=box(3)/2;
        axis(box);
    end
    ylabel('Objective');
    xlabel('Iteration');
    
    
end

function []=ExploreParamEffect_SingleStatic(optimList,optimOpts,resstruct,plotType,resDir)
    
    isField=cellfun(@isempty,regexp(fieldnames(optimOpts),'num'));
    fields=fieldnames(optimOpts);
    fields=fields(isField);
    
    for ii=1:numel(optimList)
        
        currList=[optimOpts(:).([fields{ii},'num'])];
        for jj=1:max(currList)
            figobj=struct('h',[],'axh',[]);
            figobj.h=figure;
            figobj.axh=axes;
            hold on
            colorOrder=get(gca,'colororder');
            ll=1;
            for kk=find(currList==jj)
                figobj.color=colorOrder(mod(ll-1,size(colorOrder,1))+1,:);
                figobj=PlotIndividualRefinement(resstruct(kk),plotType,figobj);
                figobj.l(1,1).DisplayName=regexprep(resstruct(kk).name{1},'_',' ');
                legEntry(ll)=figobj.l(1,1);
                ll=ll+1;
            end
            figobj.leg=legend(legEntry(1:ll-1));
            figobj.leg.Location='northeastoutside';
            if ischar(optimOpts(kk).(fields{ii}))
                figName=[fields{ii},'_',optimOpts(kk).(fields{ii})];
            elseif isnumeric(optimOpts(kk).(fields{ii})) || islogical(optimOpts(kk).(fields{ii}))
                figName=[fields{ii},'_',num2str(optimOpts(kk).(fields{ii}))];
            end
            figName=[figName,'_',plotType];
            figobj.h.Name=figName;
            hgsave(figobj.h,[resDir,filesep,'FigParam_',figName,'.fig']);
        end
        
        
    end
    
    
end

function []=ExploreParamEffect_SingleMoving(optimList,optimOpts,resstruct,plotType,resDir)
    
    isField=cellfun(@isempty,regexp(fieldnames(optimOpts),'num'));
    fields=fieldnames(optimOpts);
    fields=fields(isField);
    MFigMax=3;
    NFigMax=3;
    fullList=zeros([numel(optimOpts),numel(optimList)]);
    for ii=1:numel(optimList)
        fullList(:,ii)=[optimOpts(:).([fields{ii},'num'])]';
    end
    
    for ii=1:numel(optimList)
        actfields=fields([1:ii-1,ii+1:end]);
        currList=fullList(:,[1:ii-1,ii+1:end]);
        tightList=RemoveIdenticalVectors(currList);
        
        nFig=ceil(sqrt(size(tightList,1)));
        mFig=ceil(size(tightList,1)/nFig);
        numFig=1;
        if nFig>NFigMax && mFig>MFigMax
            nFig=NFigMax;
            mFig=MFigMax;
            numFig=size(tightList,1)/(MFigMax*NFigMax);
        end
        for figLoop=1:numFig
            figobj=struct('h',[],'axh',[]);
            figobj.h=figure;
            figList=((figLoop-1)*(MFigMax*NFigMax)+1):min(figLoop*(MFigMax*NFigMax),size(tightList,1));
            for jj=figList
                
                figobj.axh=subplot(mFig,nFig,mod(jj-1,(MFigMax*NFigMax))+1);
                hold on
                colorOrder=get(gca,'colororder');
                ll=1;
                for kk=find(all(currList==ones([size(currList,1),1])*tightList(jj,:),2))'
                    figobj.color=colorOrder(mod(ll-1,size(colorOrder,1))+1,:);
                    figobj=PlotIndividualRefinement(resstruct(kk),plotType,figobj);
                    figobj.l(1,1).DisplayName=regexprep(resstruct(kk).name{1},'_',' ');
                    if ischar(optimOpts(kk).(fields{ii}))
                        figobj.l(1,1).DisplayName=optimOpts(kk).(fields{ii});
                    elseif isnumeric(optimOpts(kk).(fields{ii})) || islogical(optimOpts(kk).(fields{ii}))
                        figobj.l(1,1).DisplayName=num2str(optimOpts(kk).(fields{ii}));
                    end
                    legEntry(ll)=figobj.l(1,1);
                    ll=ll+1;
                end
                
                figobj.leg=legend(legEntry(1:ll-1));
                if ll-1>4
                    figobj.leg.Location='northeastoutside';
                end
                nTitCell=1;
                plotName=cell(0);
                plotName{nTitCell}='Constant : ';
                for mm=1:numel(actfields)
                    if ischar(optimOpts(kk).(actfields{mm}))
                        plotName{nTitCell}=[plotName{nTitCell},optimOpts(kk).(actfields{mm}),' (',actfields{mm},')',', '];
                    elseif isnumeric(optimOpts(kk).(actfields{mm})) || islogical(optimOpts(kk).(actfields{mm}))
                        plotName{nTitCell}=[plotName{nTitCell},num2str(optimOpts(kk).(actfields{mm})),' (',actfields{mm},')',', '];
                    end
                    if numel(plotName{nTitCell})>30
                        nTitCell=nTitCell+1;
                        plotName{nTitCell}='';
                    end
                end
                title(plotName)
                
                
                
            end
            figobj.h.Name=[fields{ii},'_',plotType,'_',int2str(figLoop),'_static'];
            hgsave(figobj.h,[resDir,filesep,'FigParam_',fields{ii},'_',plotType,'_',int2str(figLoop),'.fig']);
        end
    end
    
    
end

function []=ExtractCarpetParamPair(paramPair,optimList,optimOpts,resstruct,resDir)
    
    
    isField=cellfun(@isempty,regexp(fieldnames(optimOpts),'num'));
    fields=fieldnames(optimOpts);
    fields=fields(isField);
    MFigMax=3;
    NFigMax=3;
    fullList=zeros([numel(optimOpts),numel(optimList)]);
    for ii=1:numel(optimList)
        fullList(:,ii)=[optimOpts(:).([fields{ii},'num'])]';
    end
    for ii=1:numel(paramPair)
        pairInOpts(ii)=find(cellfun(@isempty,regexprep(fields,paramPair{ii},'')));
    end
    
    trimList=fullList;trimList(:,pairInOpts)=[];
    shortOptList=RemoveIdenticalVectors(trimList);
    actfields=fields;
    actfields(pairInOpts)=[];
    for ii=1:numel(optimList)
        if isnumeric(optimList{ii}{1})
            actOpts=vertcat(optimList{ii}{:});
        elseif ischar(optimList{ii}{1})
            
            actOpts=char(optimList{ii}{:});
        else
            actOpts=[1:numel(optimList{ii})]';
        end
        [~,newOrd{ii}]=sort(actOpts(:,1));
    end
    for ii=1:size(shortOptList,1)
        actRefines=find(all(trimList==ones([size(trimList,1),1])*shortOptList(ii,:),2))';
        allMembers=fullList(actRefines,pairInOpts);
        maxRef=0;
        kk=1;
        for jj=actRefines
            actMember=[[resstruct(jj).refine(:).stage];[resstruct(jj).refine(:).objend]];
            maxRef=max(maxRef,max(actMember(1,:)));
            actMember=actMember(:);
            allMembers(kk,(numel(pairInOpts)+1):(numel(pairInOpts)+numel(actMember)))...
                =actMember;
            kk=kk+1;
        end
        figobj=struct('h',[]);
        figobj.h=figure;
        figobj.axh=axes;
        colorOrder=get(gca,'colororder');
        hold on
        for jj=0:maxRef
            figobj.color=colorOrder(mod(jj-1,size(colorOrder,1))+1,:);
            figobj=PlotCarpetPlot(allMembers(:,[1:numel(pairInOpts),...
                numel(pairInOpts)+(jj+1)*2]),optimList(pairInOpts),newOrd(pairInOpts),figobj);
            figobj.l3.DisplayName=['Ref ',int2str(jj)];
            lEntry(jj+1)=figobj.l3;
        end
        figobj.leg=legend(lEntry(1:maxRef+1));
        figobj.leg.Location='northeastoutside';
        
        
        kk=actRefines(1);
        nTitCell=1;
        plotName=cell(0);
        plotName{nTitCell}='Constant : ';
        for mm=1:numel(actfields)
            if ischar(optimOpts(kk).(actfields{mm}))
                plotName{nTitCell}=[plotName{nTitCell},optimOpts(kk).(actfields{mm}),' (',actfields{mm},')',', '];
            elseif isnumeric(optimOpts(kk).(actfields{mm})) || islogical(optimOpts(kk).(actfields{mm}))
                plotName{nTitCell}=[plotName{nTitCell},num2str(optimOpts(kk).(actfields{mm})),' (',actfields{mm},')',', '];
            end
            if numel(plotName{nTitCell})>30
                nTitCell=nTitCell+1;
                plotName{nTitCell}='';
            end
        end
        title(plotName)
        
        figobj.h.Name=[paramPair{1},'-',paramPair{2},'_',int2str(ii)];
        hgsave(figobj.h,[resDir,filesep,'CarpetPlot_',figobj.h.Name,'.fig']);
    end
    
    
    
end


function [figobj]=PlotCarpetPlot(data,labels,order,figobj)
    
    
    isLoadColor=false;
    if ~exist('figobj','var')
        figobj.h=figure('Name',resstruct.name{1});
        figobj.axh=axes;
        isLoad=false;
    else
        figure(figobj.h)
        axes(figobj.axh)
        isLoadColor=isfield(figobj,'color');
        isLoad=true;
    end
    
    x=1:max(max(data(:,1)),numel(labels{1}));
    y=1:max(max(data(:,2)),numel(labels{2}));
    
    [X,Y]=meshgrid(x,y);
    
    Z=ones(size(X))*00;
    
    for ii=1:size(data,1)
        Z(data(ii,2),data(ii,1))=data(ii,3);
    end
    Z=Z(order{2},order{1});
    figobj.mesh=mesh(X,Y,Z);
    figobj.l3=plot3(FindObjNum([],data(:,1),order{1}),FindObjNum([],data(:,2),order{2}),data(:,3));
    figobj.l3.LineStyle='none';
    figobj.l3.Marker='*';
    if isLoadColor
        figobj.mesh.CData=repmat(reshape(figobj.color,[1 1 3]),size(figobj.mesh.CData));
        figobj.l3.Color=figobj.color;
    end
    figobj.axh.ZScale='log';
    figobj.axh.XTick=1:numel(labels{1});
    figobj.axh.YTick=1:numel(labels{2});
    figobj.axh.XTickLabel=labels{1}(order{1});
    figobj.axh.XTickLabelRotation=45;
    figobj.axh.YTickLabelRotation=45;
    figobj.axh.YTickLabel=labels{2}(order{2});
    hidden off
    
    
end


%% Load Data

function [refinestruct]=TreatPathsIntoStruct(pltPaths,compareopts)
    
    disp('This code needs to be modified for function maximisation')
    
    for ii=1:numel(pltPaths)
        resstruct(ii)=TreatPathData(pltPaths{ii},compareopts);
        
    end
    resstruct=resstruct(~cellfun(@isempty,{resstruct(:).obj}));
    compList=[{'pattern'},compareopts.optim,compareopts.snake];
    kk=1;
    for ii=1:numel(compList)
        if iscell(compList{kk})
            kk=kk+numel(compList{ii})-1;
            compList=[compList(1:ii-1),compList{ii},compList(ii+1:end)];
        end
        kk=kk+1;
    end
    [resstruct,patternList]=IdentifyUniqueOptions(resstruct,compList);
    %patternList=patternList{1};
    pattNumRes=zeros([numel(compList),numel(resstruct)]);
    for ii=1:numel(compList)
        
        pattNumRes(ii,:)=[resstruct(:).([compList{ii},'num'])];
        
    end
    pattNumRes=pattNumRes';
    pattNumRes2=pattNumRes*(cumsum(max(pattNumRes))');
    
    [pattNumRes2,posPatNumRes2,allPatNumRes]=unique(pattNumRes2);
    
    for ii=numel(pattNumRes2):-1:1
        for jj=1:size(pattNumRes,2)
            refinestruct(ii).name{jj}=patternList{jj}{pattNumRes(posPatNumRes2(ii),jj)};
        end
        refinestruct(ii).refine=resstruct(pattNumRes2(allPatNumRes)==pattNumRes2(ii));
        [~,sortStruct]=sort([refinestruct(ii).refine(:).stage]);
        refinestruct(ii).refine=refinestruct(ii).refine(sortStruct);
        [refinestruct(ii).refine]=TrimRefineStruct(refinestruct(ii).refine);
        refinestruct(ii).paramoptim=refinestruct(ii).refine(1).paramoptim;
        refinestruct(ii).paramsnake=refinestruct(ii).refine(1).paramsnake;
        dirPath=regexp(refinestruct(ii).refine(1).path,filesep);
        refinestruct(ii).path=refinestruct(ii).refine(1).path(1:dirPath(end-1));
    end
    
end

function [refinestruct]=TreatPathsIntoStructOld(pltPaths,compareopts)
    
    disp('This code needs to be modified for function maximisation')
    
    for ii=1:numel(pltPaths)
        resstruct(ii)=TreatPathData(pltPaths{ii},compareopts);
        
    end
    resstruct=resstruct(~cellfun(@isempty,{resstruct(:).obj}));
    [resstruct,patternList]=IdentifyUniqueOptions(resstruct,{'pattern'});
    patternList=patternList{1};
    pattNumRes=[resstruct(:).patternnum];
    for ii=numel(patternList):-1:1
        refinestruct(ii).name=patternList{ii};
        refinestruct(ii).refine=resstruct(pattNumRes==ii);
        [~,sortStruct]=sort([refinestruct(ii).refine(:).stage]);
        refinestruct(ii).refine=refinestruct(ii).refine(sortStruct);
        [refinestruct(ii).refine]=TrimRefineStruct(refinestruct(ii).refine);
        refinestruct(ii).paramoptim=refinestruct(ii).refine(1).paramoptim;
        refinestruct(ii).paramsnake=refinestruct(ii).refine(1).paramsnake;
        dirPath=regexp(refinestruct(ii).refine(1).path,filesep);
        refinestruct(ii).path=refinestruct(ii).refine(1).path(1:dirPath(end-1));
    end
    
end

function [refine]=TrimRefineStruct(refine)
    
    for ii=2:numel(refine)
        refine(ii).iterstart=refine(ii-1).iterend;
    end
    for ii=1:numel(refine)
        refine(ii).obj=refine(ii).obj(refine(ii).iterstart:refine(ii).iterend);
        refine(ii).des=refine(ii).des(refine(ii).iterstart:refine(ii).iterend,:);
        refine(ii).add=refine(ii).add(refine(ii).iterstart:refine(ii).iterend);
        refine(ii).objstart=refine(ii).obj(1);
        refine(ii).objend=refine(ii).obj(end);
    end
    
    
end

function [refinestruct]=TreatPathData(optimPath,compareopts)
    
    optimCell=regexp(optimPath,filesep,'split');
    
    dirPath=regexp(optimPath,filesep);
    dirPath=optimPath(1:dirPath(end));
    [paramPath]=FindFinalParamFiles({dirPath});
    
%     refinestruct.stage=str2double(regexprep(optimCell{end-1},['Dir.*',optimCell{end-2},'_'],''));
%     refinestruct.stage=str2double(regexp(optimCell{end-1},['_[0-9]+',filesep],'match'));
    stageRef=regexp(optimCell{end-1},'_','split');
    refinestruct.stage=str2double(stageRef{end});
    refinestruct.time=datenum(regexprep(regexprep(optimCell{end-1},'Dir_',''),...
        [optimCell{end-2},'_.*'],''),'yyyy-mm-ddTHHMMSS');
    refinestruct.pattern=optimCell{end-2};
    
    
    instruct=load(optimPath);
    [refinestruct.iterstart,refinestruct.iterend,refinestruct.ndesvar,...
        refinestruct.obj,refinestruct.des,refinestruct.add]=TrimDataToStruct(instruct.optimstruct);
    instruct=load(paramPath{1});
    
    for ii=1:numel(compareopts.optim)
        refinestruct.paramoptim.(compareopts.optim{ii})=ExtractVariables(...
            compareopts.optim(ii),instruct.paramoptim);
        refinestruct.(compareopts.optim{ii})=ExtractVariables(...
            compareopts.optim(ii),instruct.paramoptim);
        if iscell(refinestruct.paramoptim.(compareopts.optim{ii}))
           refinestruct.paramoptim.(compareopts.optim{ii})=refinestruct.paramoptim.(compareopts.optim{ii}){1};
           refinestruct.(compareopts.optim{ii})=refinestruct.paramoptim.(compareopts.optim{ii}); 
        end
    end
    if numel(compareopts.optim)==0
        refinestruct.paramoptim=struct([]);
    end
    for ii=1:numel(compareopts.snake)
        refinestruct.paramsnake.(compareopts.snake{ii})=ExtractVariables(...
            compareopts.snake(ii),instruct.paramoptim.parametrisation);
        refinestruct.(compareopts.snake{ii})=ExtractVariables(...
            compareopts.snake(ii),instruct.paramoptim.parametrisation);
        
        if iscell(refinestruct.paramsnake.(compareopts.snake{ii}))
           refinestruct.paramoptim.(compareopts.snake{ii})=refinestruct.paramoptim.(compareopts.snake{ii}){1};
           refinestruct.(compareopts.snake{ii})=refinestruct.paramoptim.(compareopts.snake{ii});
        end
    end
    if numel(compareopts.snake)==0
        refinestruct.paramsnake=struct([]);
    end
    
    refinestruct.objstart=refinestruct.obj(refinestruct.iterstart);
    refinestruct.objend=refinestruct.obj(refinestruct.iterend);
    refinestruct.path=dirPath;
    refinestruct.dir=optimCell{end-1};
end

function [iterstart,iterend,nDesVar,objHist,desHist,addHist]=TrimDataToStruct(optimstruct)
    
    iterstart=1;
    n=1;
    while n<=numel(optimstruct) && ~isempty(optimstruct(n).population(1).objective)
        n=n+1;
    end
    iterend=n-1;
    nDesVar=numel(optimstruct(iterend).population(1).fill);
    for ii=iterend:-1:1
        
        [objHist(ii),currOpt]=min([optimstruct(ii).population(:).objective]);
        desHist(ii,1:numel(optimstruct(ii).population(currOpt).fill))=...
            optimstruct(ii).population(currOpt).fill;
        addHist(ii)=optimstruct(ii).population(currOpt).additional;
    end
    
    
end

function [refinestruct,patternList]=IdentifyUniqueOptions(refinestruct,jobfields)
    
    if nargin<2
        jobfields=fieldnames(refinestruct);
    end
    
    %
    
    n=numel(refinestruct);
    
    patternList=cell(0);
    for ii=1:numel(jobfields);
        patternList{ii}={};
        
        for jj=1:n
            optfound=false;
            
            for kk=1:numel(patternList{ii})
                test=false;
                if ischar(refinestruct(jj).(jobfields{ii}))
                    test=strcmp(patternList{ii}{kk},refinestruct(jj).(jobfields{ii}));
                elseif isnumeric(refinestruct(jj).(jobfields{ii})) || islogical(refinestruct(jj).(jobfields{ii}))
                    test=(patternList{ii}{kk}==refinestruct(jj).(jobfields{ii}));
                end
                if test
                    refinestruct(jj).([jobfields{ii},'num'])=kk;
                    optfound=true;
                    break
                end
            end
            if ~optfound
                patternList{ii}{end+1}=refinestruct(jj).(jobfields{ii});
                refinestruct(jj).([jobfields{ii},'num'])=numel(patternList{ii});
            end
        end
    end
    
    
    
    
end

function [allOpts,allOptsList]=ConcatenateStructures(optimOpts,optimList,snakeOpts,snakeList)
    
    if ~isempty(optimOpts) && ~isempty(snakeOpts)
        allOpts=catstruct(optimOpts,snakeOpts);
        allOptsList=[optimList,snakeList];
    elseif ~isempty(optimOpts)
        allOpts=optimOpts;
        allOptsList=optimList;
    elseif ~isempty(snakeOpts)
        allOpts=snakeOpts;
        allOptsList=snakeList;
    else
        allOpts=struct([]);
        allOptsList=cell(0);
    end
    
    
end
%% Explore folders

function [pltPaths]=FindOptimResFiles(childFolder)
    kk=1;
    for ii=1:length(childFolder)
        intermPath=FindDir(childFolder{ii},'OptimRes',false);
        if ~isempty(intermPath)
            if numel(intermPath)>1 && ~isempty(find(cellfun(@isempty,regexp(intermPath,...
                    'partial')),1,'first'))
                intermPath=intermPath(find(cellfun(@isempty,regexp(intermPath,...
                    'partial')),1,'first'));
            end
            pltPaths(kk)=intermPath(1);
            kk=kk+1;
        end
    end
    
end

function [pltPaths]=FindFinalParamFiles(childFolder)
    kk=1;
    for ii=1:length(childFolder)
        intermPath=FindDir(childFolder{ii},'FinalParam',false);
        if ~isempty(intermPath)
            
            if numel(intermPath)>1 && ~isempty(find(cellfun(@isempty,regexp(intermPath,...
                    'partial')),1,'first'))
                intermPath=intermPath(find(cellfun(@isempty,regexp(intermPath,...
                    'partial')),1,'first'));
            end
            pltPaths(kk)=intermPath(1);
            kk=kk+1;
        end
    end
    
end

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    if ~isempty(subDir)
        for ii=1:length(subDir)
            subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
                ~xor(subDir(ii).isdir,isTargDir);
        end
        
        returnSub=find([subDir(:).isProfile]);
    else
        returnSub=[];
    end
    
    if isempty(returnSub)
        %disp('FindDir Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end

function [addFolders]=ExploreFolderTree(rootDir)
    % adds a set of paths to the active path
    addFolders=rootDir;
    for ii=1:length(rootDir)
        dirinfo=dir(rootDir{ii});
        dirNames={dirinfo([dirinfo(:).isdir]).name};
        dirNames(1:2)=[];
        if numel(dirNames)>0
            branchDir={''};
            for jj=1:length(dirNames)
                branchDir{jj}=[rootDir{ii},filesep,dirNames{jj}];
            end
            
            [addSubFolders]=ExploreFolderTree(branchDir);
            addFolders=[addFolders,addSubFolders];
        end
    end
    
end

