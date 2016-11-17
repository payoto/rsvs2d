function [dirDat,sortDat]=CompareAllGradients(rootDir,maxRecurs)
    
    [childFolder]=ExploreFolderTree(rootDir,maxRecurs);
    [resDir]=FindResDirFolder(childFolder);
    
    
    [dirDat]=PostTreatFiles(resDir);
    
    dirSummarryFigs=[rootDir{1},filesep,'GradientComp'];
    dateStr=datestr(now,'yy_mm_dd_HH_MM');
    mkdir(dirSummarryFigs)
    save([dirSummarryFigs,filesep,'AllParam_',dateStr,'.mat'],'dirDat')
    [sortDat]=AnalyseGradients(dirDat,dirSummarryFigs);
end

function [h,directionChange]=PlotGradients(paramoptim,optimstruct)
    
    normVec=@(vec) sqrt(sum(vec.^2,2));
    normaliseArray=@(array)  array./(sqrt(sum(array.^2,2))*ones(1,size(array,2))); 
    
    supportOptim=paramoptim.optim.supportOptim;
    
    h=figure('Name',['Gradients_',ExtractVariables({'optimCase'},paramoptim)]...
        ,'Position',[ 100 150 1400 700]);
    symDesVarList=ExtractVariables({'symDesVarList'},paramoptim);
    rmCol=symDesVarList(2,:);
    subplot(2,2,1)
    grads=vertcat(supportOptim.hist(:).gradfk);
    grads(:,rmCol)=[];
    surf(log10(abs(grads)))
    ylabel('iteration')
    xlabel('design variable')
    title('gradients')
    
    view(0,90)
    subplot(2,2,2)
    grads=vertcat(supportOptim.hist(:).prevDir);
    grads(:,rmCol)=[];
    surf(log10(abs(grads)))
    ylabel('iteration')
    xlabel('design variable')
    title('previous direction')
    view(0,90)
    subplot(2,2,3)
    grads=vertcat(supportOptim.hist(:).prevStep);
    surf(((grads)),(abs((grads))))
    ylabel('iteration')
    xlabel('design variable')
    title('previous direction')
    view(0,90)
    subplot(2,2,4)
    hold on
    grads=vertcat(supportOptim.hist(:).gradfk);
    grads(:,rmCol)=[];
    gradsm1=vertcat(supportOptim.hist(:).gradfkm1);
    gradsm1(:,rmCol)=[];
    directionChange.Grad=dot(normaliseArray(grads),normaliseArray(gradsm1),2);
    directionChange.Grad(1)=1;
    gradStep=vertcat(supportOptim.hist(:).prevStep);
    gradStep(:,rmCol)=[];
    directionChange.Step=dot(normaliseArray(grads),normaliseArray(gradStep),2);
    gradDir=vertcat(supportOptim.hist(:).prevDir);
    gradDir(:,rmCol)=[];
    directionChange.Dir=dot(normaliseArray(grads),normaliseArray(gradDir),2);
    directionChange.grads=grads;
    directionChange.gradsm1=gradsm1;
    kk=1;
    fillPrec=zeros(size(optimstruct(1).population(1).fill));
    for ii=1:2:numel(optimstruct)
        changePos(kk)=normVec(fillPrec-optimstruct(ii).population(1).fill);
        fillInf(kk,1:numel(optimstruct(ii).population(1).fill))...
            =optimstruct(ii).population(1).fill;
        fillPrec=optimstruct(ii).population(1).fill;
        kk=kk+1;
    end
    changePos(1)=ExtractVariables({'startVol'},paramoptim);
    directionChange.Pos=changePos;
    
    step=1:numel(supportOptim.hist);
    ii=1;
    l(ii)=plot(step,directionChange.Grad);
    l(ii).DisplayName='Change in direction - gradient';
    ii=ii+1;
    l(ii)=plot(step,directionChange.Step);
    l(ii).DisplayName='gradient // direction';
    ii=ii+1;
    l(ii)=plot(step,directionChange.Dir);
    l(ii).DisplayName='gradient // step';
    ii=ii+1;
    try
        scaleEvol=[supportOptim.hist(:).scale];
        l(ii)=plot(step,scaleEvol);
        l(ii).DisplayName='scale';
        ii=ii+1;
    catch
    end
    try
        iterEvol=[supportOptim.hist(:).iter];
        l(ii)=plot(step,iterEvol);
        l(ii).DisplayName='iter since last refresh';
        ii=ii+1;
    catch
    end
    l(ii)=plot(step,changePos);
    l(ii).DisplayName='Length of movement';
    ii=ii+1;
    legend(l)
    
end


%% Exploration functions


function [resDir]=FindResDirFolder(childFolder)
    
    isResDir=false(size(childFolder));
    
    for ii=1:numel(childFolder)
        foldName=regexp(childFolder{ii},filesep,'split');
        isResDir(ii)=~isempty(regexp(foldName{end},'Dir_', 'once'));
        
    end
    resDir=childFolder(isResDir);
end

function [postDir]=NeedsPostTreat(resDir)
    
    isPost=false(size(resDir));
    
    for ii=1:numel(resDir)
        
        isPost(ii)=isempty(FindDir(resDir{ii},'Optimal__',true)) || ...
            isempty(FindDir(resDir{ii},'GradHistory_',false));
        
    end
    postDir=resDir(isPost);
end

function [dirDat]=PostTreatFiles(postDir)
    
    
    fprintf('\n%i folders to be postreated \n',numel(postDir))
    for ii=1:numel(postDir)
        fprintf('%s\n',postDir{ii})
    end
    fprintf('\n')
    dirDat=repmat(struct('path','','case','','supportOptim',struct([]),...
        'diffStep',[],'paramoptim',struct([])),[1,numel(postDir)]);
    rmII=[];
    for ii=1:numel(postDir)
        
        dirDat(ii).path=postDir{ii};
        
        try 
            [dirDat(ii).case,dirDat(ii).supportOptim,dirDat(ii).diffStep,dirDat(ii).paramoptim]=...
            ExtractParamOptim(postDir{ii});
            [optimstruct]=ExtractOptimStruct(postDir{ii});
            [h,dirDat(ii).dirChange]=PlotGradients(dirDat(ii).paramoptim,optimstruct);
            hgsave(h,[postDir{ii},filesep,'GradHistory_',h.Name])
            close(h)
        catch
            rmII=[rmII,ii];
        end
    end
    
    if ~isempty(rmII)
        dirDat(rmII)=[];
        warning([int2str(numel(rmII)),' Failed:'])
        char(postDir{rmII})
    end
end

function [T]=CaptureOutPostTreatIncomplete(postDir,optimstruct)
    [T,~]=evalc('PostTreatIncomplete(postDir,[],optimstruct)');
end

function [optimstruct]=ExtractOptimStruct(postDir)
    
    [MatPath,MatName]=FindDir(postDir,'OptimRes',false);
    if isempty(regexp(MatName{1},'partial','once'))
        load(MatPath{1},'optimstruct')
    else
        load(MatPath{1},'optimstruct')
        
    end
    maxIter=numel(optimstruct);
    while isempty(optimstruct(maxIter).population(1).objective)
        maxIter=maxIter-1;
    end
    optimstruct=optimstruct(1:maxIter);
end

function [caseName,supportOptim,diffStep,paramoptim]=ExtractParamOptim(postDir)
    
    [MatPath,MatName]=FindDir(postDir,'FinalParam',false);
    if ~isempty(MatName)
        if isempty(regexp(MatName{1},'partial','once'))
            load(MatPath{1})
        else
            load(MatPath{1})

        end
    end
    
    if exist('paramoptim','var')
        supportOptim=paramoptim.optim.supportOptim;
        caseName=paramoptim.general.optimCase;
        diffStep=paramoptim.optim.CG.minDiffStep;
    else
        error('paramoptim was not loaded from the .mat file')
    end
    
    %     maxIter=numel(optimstruct);
    %     while isempty(optimstruct(maxIter).population(1).objective)
    %         maxIter=maxIter-1;
    %     end
    %     optimstruct=optimstruct(1:maxIter);
end

%

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
        disp('FindDir Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end

function [addFolders]=ExploreFolderTree(rootDir,maxRecurs)
    % adds a set of paths to the active path
    if nargin<2;maxRecurs=inf;end
    
    addFolders=rootDir;
    maxRecurs=maxRecurs-1;
    for ii=1:length(rootDir)
        dirinfo=dir(rootDir{ii});
        dirNames={dirinfo([dirinfo(:).isdir]).name};
        dirNames(1:2)=[];
        if numel(dirNames)>0
            branchDir={''};
            for jj=1:length(dirNames)
                branchDir{jj}=[rootDir{ii},filesep,dirNames{jj}];
            end
            if maxRecurs>0
                
                [addSubFolders]=ExploreFolderTree(branchDir,maxRecurs);
                addFolders=[addFolders,addSubFolders];
            end
        end
    end
    
end