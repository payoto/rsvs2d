function [resstruct]=PostTreatGroupRun(entryPoint,varargin)
    % This function post treats an entire group of runs.
    
    switch entryPoint
        case 'data'
            resstruct=CollectGroupRunInformation(varargin{:});
        case 'plot'
            resstruct=varargin{1};
    end
    
    PlotGroupedInformation(resstruct)
end

function [resstruct]=CollectGroupRunInformation(rootDir,compareopts)
    
    [childFolder]=ExploreFolderTree({rootDir});
    [optimPaths]=FindOptimResFiles(childFolder);
    %[paramPaths]=FindFinalParamFiles(childFolder);
    
    [resstruct]=TreatPathsIntoStruct(optimPaths,compareopts);
end

function []=PlotGroupedInformation(resstruct)
    for ii=1:numel(resstruct)
        PlotIndividualRefinement(resstruct(ii))
    end
    
end

%% Plot Data

function []=PlotIndividualRefinement(resstruct,lcol,lpat)
    
    
    figobj.h=figure('Name',resstruct.name);
    figobj.axh=axes;
    hold on
    for ii=1:numel(resstruct.refine)
        
        figobj.l(ii,1)=plot(resstruct.refine(ii).iterstart:resstruct.refine(ii).iterend,...
            resstruct.refine(ii).obj);
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
    figobj.axh.YScale='log';
    box=axis;
    box(2)=box(2)+10;
    %box(3)=box(3)/2;
    axis(box);
    ylabel('Objective');
    xlabel('Iteration');
    
    hgsave(figobj.h,[resstruct.path,filesep,'FigRefOpt_',resstruct.name,'.fig']);
end


%% Load Data

function [refinestruct]=TreatPathsIntoStruct(pltPaths,compareopts)
    
    disp('This code needs to be modified for function maximisation')
    
    for ii=1:numel(pltPaths)
        resstruct(ii)=TreatPathData(pltPaths{ii},compareopts);
        
    end
    [resstruct,patternList]=IdentifyUniqueOptions(resstruct);
    
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
    iterend=numel(optimstruct);
    nDesVar=numel(optimstruct(end).population(1).fill);
    
    for ii=numel(optimstruct):-1:1
        
        [objHist(ii),currOpt]=min([optimstruct(ii).population(:).objective]);
        desHist(ii,1:numel(optimstruct(ii).population(currOpt).fill))=...
            optimstruct(ii).population(currOpt).fill;
        addHist(ii)=optimstruct(ii).population(currOpt).additional;
    end
    
    
end

function [refinestruct,patternList]=IdentifyUniqueOptions(refinestruct)
    
    %jobfields=fieldnames(refinestruct);
    jobfields={'pattern'};
    n=numel(refinestruct);
    
    
    ii=1;
    patternList={};
    
    for jj=1:n
        optfound=false;
        
        for kk=1:numel(patternList)
            if strcmp(patternList{kk},refinestruct(jj).(jobfields{ii}))
                refinestruct(jj).patternnum=kk;
                optfound=true;
                break
            end
        end
        if ~optfound
            patternList{end+1}=refinestruct(jj).(jobfields{ii});
            refinestruct(jj).patternnum=numel(patternList);
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

