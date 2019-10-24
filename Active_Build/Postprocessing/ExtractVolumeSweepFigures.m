function [resstruct]=ExtractVolumeSweepFigures(pathStr)

    [childFolder]=ExploreFolderTree({pathStr});
    [optimPaths]=FindOptimResFiles(childFolder);
    %[paramPaths]=FindFinalParamFiles(childFolder);
    compareOpts=struct('optim',{{'desVarVal',{'resVal',2}}},'snake',{{'axisRatio'}});
    [resstruct]=TreatPathsIntoStruct(optimPaths,compareOpts);

    plotstruct=ReturnLineProfile(resstruct);
    [~,sOrd]=sort([plotstruct(:).A]);
    plotstruct=plotstruct(sOrd);
    
    [h(1)]=PlotLine(plotstruct);
    h(2)=PlotProfiles(plotstruct);
    h(3)=PlotProfilesScaled(plotstruct);
    SaveFigs(h,pathStr)
end

function [plotstruct]=ReturnLineProfile(resstruct)
    for ii=1:numel(resstruct)
        flag=true;
        plotstruct(ii).A=resstruct(ii).paramoptim.desVarVal;
        plotstruct(ii).axisRatio=resstruct(ii).paramsnake.axisRatio;
        while flag && numel(resstruct(ii).refine)>0
            [plotstruct(ii).obj,loc]=min([resstruct(ii).refine(:).objend]);
            loopPath=FindDir(resstruct(ii).refine(loc).optimSolPath,'restart',0);
            if numel(loopPath)>0
                instruct=load(loopPath{1},'loop');
                for jj=1:numel(instruct.loop)
                    plotstruct(ii).points{jj}=instruct.loop(jj).subdivision;
                end
                flag=false;
            else
                resstruct(ii).refine(loc)=[];
                flag=true;
            end
        end
        
    end
    
    
end

function [dy]=AdaptSizeforBusemann(e)
    include_Optimisation
    [loop]=ConstantArea_Busemann(0,1,e,2);
    ymax=nan;
    ymin=nan;
    for ii=1:numel(loop)
        ymax=max([loop(ii).subdivision(:,2);ymax]);
        ymin=min([loop(ii).subdivision(:,2);ymin]);
    end
    dy=(ymax-ymin)/2;
end

function h=PlotProfiles(plotstruct)
    
    cOrd=get(gca,'ColorOrder');
    lStyle={'-','--'};
    
    h=figure('Name','Profile Group');
    hold on
    maxT=zeros(numel(plotstruct),2);
    for ii=1:numel(plotstruct)
        
        for jj=1:numel(plotstruct(ii).points)
            l(ii)=plot(plotstruct(ii).points{jj}(:,1),plotstruct(ii).points{jj}(:,2), ...
                'LineStyle',lStyle{mod(floor((ii-1)/size(cOrd,1)),2)+1},...
                'Color',cOrd(mod(ii-1,size(cOrd,1))+1,:),...
                'DisplayName',num2str(plotstruct(ii).A,'%.2f'));
            if numel(plotstruct(ii).points)==1
                [~,posMaxT]=max(plotstruct(ii).points{jj}(:,2));
                
                
                plot(plotstruct(ii).points{jj}(posMaxT,1),plotstruct(ii).points{jj}(posMaxT,2), ...
                    'LineStyle',lStyle{mod(floor((ii-1)/size(cOrd,1)),2)+1},...
                    'Color',cOrd(mod(ii-1,size(cOrd,1))+1,:),...
                    'DisplayName',num2str(plotstruct(ii).A,'%.2e'),...
                    'Marker','*');
                
                maxT(ii,:)=[plotstruct(ii).points{jj}(posMaxT,1),plotstruct(ii).points{jj}(posMaxT,2)];
            end
        end
    end
    plot(maxT(:,1),maxT(:,2),'k-')
    legend(l,'Location','EastOutside')
    
end

function h=PlotProfilesScaled(plotstruct)
    
    cOrd=get(gca,'ColorOrder');
    lStyle={'-','--'};
    
    h=figure('Name','Profile Group');
    hold on
    maxT=zeros(numel(plotstruct),2);
    for ii=1:numel(plotstruct)
        [dy]=AdaptSizeforBusemann(plotstruct(ii).A);
        for jj=1:numel(plotstruct(ii).points)
            l(ii)=plot(plotstruct(ii).points{jj}(:,1),plotstruct(ii).points{jj}(:,2)/dy, ...
                'LineStyle',lStyle{mod(floor((ii-1)/size(cOrd,1)),2)+1},...
                'Color',cOrd(mod(ii-1,size(cOrd,1))+1,:),...
                'DisplayName',num2str(plotstruct(ii).A,'%.2f'));
            if numel(plotstruct(ii).points)==1
                [~,posMaxT]=max(plotstruct(ii).points{jj}(:,2));
                
                
                plot(plotstruct(ii).points{jj}(posMaxT,1),plotstruct(ii).points{jj}(posMaxT,2)/dy, ...
                    'LineStyle',lStyle{mod(floor((ii-1)/size(cOrd,1)),2)+1},...
                    'Color',cOrd(mod(ii-1,size(cOrd,1))+1,:),...
                    'DisplayName',num2str(plotstruct(ii).A,'%.2e'),...
                    'Marker','*');
                
                maxT(ii,:)=[plotstruct(ii).points{jj}(posMaxT,1),plotstruct(ii).points{jj}(posMaxT,2)/dy];
            end
        end
    end
    plot(maxT(:,1),maxT(:,2),'k-')
    legend(l,'Location','EastOutside')
    
end

function [h]=PlotLine(plotstruct)
    
    A=[plotstruct(:).A];
    obj=[plotstruct(:).obj];
    h=figure('Name','Drag Evol');
    semilogy(A,obj)
    
end

function SaveFigs(h,rootDir)
    
    
    for ii=1:numel(h)
        
        hgsave(h,[rootDir,filesep,'Fig_',int2str(ii),'_',h(ii).Name,'.fig'])
        
    end
    
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
    for ii=numel(compList):-1:1
        if isnumeric(compList{ii})
            compList(ii) = [];
        end
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

function [refine]=TrimRefineStruct(refine)
    
    %     for ii=2:numel(refine)
    %         refine(ii).iterstart=refine(ii-1).iterend;
    %     end
    for ii=1:numel(refine)
        refine(ii).stage(isnan(refine(ii).stage))=0;
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
    
    rootPath=regexp(optimPath,filesep);
    
    %     refinestruct.stage=str2double(regexprep(optimCell{end-1},['Dir.*',optimCell{end-2},'_'],''));
    %     refinestruct.stage=str2double(regexp(optimCell{end-1},['_[0-9]+',filesep],'match'));
    stageRef=regexp(optimCell{end-1},'_','split');
    refinestruct.stage=str2double(stageRef{end});
    refinestruct.time=datenum(regexprep(regexprep(optimCell{end-1},'Dir_',''),...
        [optimCell{end-2},'_.*'],''),'yyyy-mm-ddTHHMMSS');
    refinestruct.pattern=optimCell{end-2};
    
    
    instruct=load(optimPath);
    [refinestruct.iterstart,refinestruct.iterend,refinestruct.ndesvar,...
        refinestruct.obj,refinestruct.des,refinestruct.add,...
        refinestruct.optimSolPath]=TrimDataToStruct(...
        RewriteOptimForPC(instruct.optimstruct,dirPath(1:rootPath(end-1))));
    optimSolPath=FindDir(dirPath,'Optimal',1);
    if numel(optimSolPath)>0
        refinestruct.optimSolPath=optimSolPath{1};
    end
    instruct=load(paramPath{1});
    
    fields = {'optim','snake'};
    param = {instruct.paramoptim, instruct.paramoptim.parametrisation};

    for jj = 1:2
        fieldAct = ['param',fields{jj}];
        for ii=1:numel(compareopts.(fields{jj}))
            if iscell(compareopts.(fields{jj}){ii})
                compareOpt = compareopts.(fields{jj}){ii}{1};
                compareOptNum = compareopts.(fields{jj}){ii}{2};
            else
                compareOpt = compareopts.(fields{jj}){ii};
                compareOptNum = 1;
            end
            refinestruct.(fieldAct).(compareOpt)=ExtractVariables(...
                {compareOpt},param{jj});
            refinestruct.(compareOpt)=refinestruct.(fieldAct).(compareOpt);
            if iscell(refinestruct.(fieldAct).(compareOpt))
                refinestruct.(fieldAct).(compareOpt) =...
                    refinestruct.(fieldAct).(compareOpt){compareOptNum};
                refinestruct.(compareOpt) =...
                    refinestruct.(fieldAct).(compareOpt);
            end
        end
        if numel(compareopts.(fields{jj}))==0
            refinestruct.(fieldAct)=struct([]);
        end
    end

    
    refinestruct.objstart=refinestruct.obj(refinestruct.iterstart);
    refinestruct.objend=refinestruct.obj(refinestruct.iterend-1);
    refinestruct.path=dirPath;
    refinestruct.dir=optimCell{end-1};
end

function [iterstart,iterend,nDesVar,objHist,desHist,addHist,...
        optimSolPath]=TrimDataToStruct(optimstruct)
    
    iterstart=1;
    n=1;
    while n<=numel(optimstruct) && ~isempty(optimstruct(n).population(1).objective)
        n=n+1;
    end
    iterend=n-1;
    nDesVar=numel(optimstruct(iterend).population(1).fill);
    for ii=iterend:-1:1
        
        [objHist(ii),currOpt]=min([optimstruct(ii).population(:).objective]);
        if ii==iterend
            optimSolPath=optimstruct(ii).population(currOpt).location;
        end
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
                elseif isnumeric(refinestruct(jj).(jobfields{ii})) ...
                        || islogical(refinestruct(jj).(jobfields{ii}))
                    test=(patternList{ii}{kk}==refinestruct(jj).(jobfields{ii}));
                end
                if test
                    refinestruct(jj).([jobfields{ii},'num'])=kk;
                    optfound=true;
                    break
                end
            end
            if ~optfound
                try
                    patternList{ii}{end+1}=refinestruct(jj).(jobfields{ii});
                    refinestruct(jj).([jobfields{ii},'num'])=numel(patternList{ii});
                catch
                end
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
    if numel(subDir)>0
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
            %fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
        end
        for ii=1:length(returnSub)
            returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
            returnName{ii}=subDir(returnSub(ii)).name;
            
        end
    else
        warning('Invalid Directory')
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

