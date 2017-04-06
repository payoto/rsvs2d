function [resstruct]=PostTreatGroupRun(entryPoint,varargin)
    % This function post treats an entire group of runs.
    plotStart=1;
    switch entryPoint
        case 'data'
            resstruct=CollectGroupRunInformation(varargin{:});
        case 'plot'
            resstruct=varargin{1};
            if numel(varargin)>1
                plotStart=varargin{2};
            end
    end
    
    PlotGroupedInformation(resstruct,plotStart)
end

function [resstruct]=CollectGroupRunInformation(rootDir,compareopts)
    
    [childFolder]=ExploreFolderTree({rootDir});
    [optimPaths]=FindOptimResFiles(childFolder);
    %[paramPaths]=FindFinalParamFiles(childFolder);
    
    [resstruct]=TreatPathsIntoStruct(optimPaths,compareopts);
    save([rootDir,filesep,'AllRefineResults.mat'],'resstruct')
end

function []=PlotGroupedInformation(resstruct,plotStart)
    
    
    
    if plotStart<2
        for ii=1:numel(resstruct)
            figobj=PlotIndividualRefinement(resstruct(ii),'full');
            hgsave(figobj.h,[resstruct(ii).path,filesep,'FigRefOpt_',...
                resstruct(ii).name,'.fig']);
            close(figobj.h)
        end
    end
    optimOpts=[resstruct(:).paramoptim];
    snakeOpts=[resstruct(:).paramsnake];
    [optimOpts,optimList]=IdentifyUniqueOptions(optimOpts);
    [snakeOpts,snakeList]=IdentifyUniqueOptions(snakeOpts);
    
    resDir=resstruct(1).path;
    n=regexp(resDir,filesep);
    resDir=[resDir(1:n(end-1)),'results'];
    mkdir(resDir)
    mkdir([resDir,filesep,'static'])
    
    plotType='ref';
    if plotStart<3
        ExploreParamEffect_SingleStatic(optimList,optimOpts,resstruct,plotType,[resDir,filesep,'static'])
    end
    mkdir([resDir,filesep,'moving'])
    if plotStart<4
        ExploreParamEffect_SingleMoving(optimList,optimOpts,resstruct,plotType,[resDir,filesep,'moving'])
    end
    
    if plotStart<5
        
    end
end

%% Plot Data

function [figobj]=PlotIndividualRefinement(resstruct,typePlot,figobj)
    
    isLoadColor=false;
    if ~exist('figobj','var')
        figobj.h=figure('Name',resstruct.name);
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
                figobj.l(1,1).DisplayName=regexprep(resstruct(kk).name,'_',' ');
                legEntry(ll)=figobj.l(1,1);
                ll=ll+1;
            end
            figobj.leg=legend(legEntry(1:ll-1));
            figobj.leg.Location='northeastoutside';
            if ischar(optimOpts(kk).(fields{ii}))
                figName=[fields{ii},'_',optimOpts(kk).(fields{ii})];
            elseif isnumeric(optimOpts(kk).(fields{ii}))
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
                    figobj.l(1,1).DisplayName=regexprep(resstruct(kk).name,'_',' ');
                    if ischar(optimOpts(kk).(fields{ii}))
                        figobj.l(1,1).DisplayName=optimOpts(kk).(fields{ii});
                    elseif isnumeric(optimOpts(kk).(fields{ii}))
                        figobj.l(1,1).DisplayName=num2str(optimOpts(kk).(fields{ii}));
                    end
                    legEntry(ll)=figobj.l(1,1);
                    ll=ll+1;
                end
                
                figobj.leg=legend(legEntry(1:ll-1));
                if ll-1>4
                    figobj.leg.Location='northeastoutside';
                end
                plotName='Constant Param: ';
                for mm=1:numel(actfields)
                    if ischar(optimOpts(kk).(actfields{mm}))
                        plotName=[plotName,optimOpts(kk).(actfields{mm}),', '];
                    elseif isnumeric(optimOpts(kk).(actfields{mm}))
                        plotName=[plotName,num2str(optimOpts(kk).(actfields{mm})),', '];
                    end
                end
                title(plotName)
                
                
                
            end
            figobj.h.Name=[fields{ii},'_',plotType,'_',int2str(figLoop),'_static'];
            hgsave(figobj.h,[resDir,filesep,'FigParam_',fields{ii},'_',plotType,'_',int2str(figLoop),'.fig']);
        end
    end
    
    
end

function []=PlotCarpetParamPair()
    
    
    
end
%% Load Data

function [refinestruct]=TreatPathsIntoStruct(pltPaths,compareopts)
    
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
    
    refinestruct.stage=str2double(regexprep(optimCell{end-1},['Dir.*',optimCell{end-2},'_'],''));
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
    end
    if numel(compareopts.optim)==0
        refinestruct.paramoptim=struct([]);
    end
    refinestruct.paramsnake=struct([]);
    for ii=1:numel(compareopts.snake)
        refinestruct.paramsnake.(compareopts.snake{ii})=ExtractVariables(...
            compareopts.snake(ii),instruct.paramoptim);
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
                elseif isnumeric(refinestruct(jj).(jobfields{ii}))
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

