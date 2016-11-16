function []=CompareAllGradients(rootDir,maxRecurs)
    
    [childFolder]=ExploreFolderTree(rootDir,maxRecurs);
    [resDir]=FindResDirFolder(childFolder);
    
    
    [dirDat]=PostTreatFiles(postDir);
    
    
end




function [h]=PlotGradients(paramoptim)
    
    supportOptim=paramoptim.optim.supportOptim;
    
    h=figure('Name',['Gradients_',ExtractVariables({'optimCase'},paramoptim(ii))]...
        ,'Position',[ 100 150 1400 700]);
    subplot(2,2,1)
    grads=vertcat(supportOptim.hist(:).gradfk);
    surf(log10(abs(grads)))
    ylabel('iteration')
    xlabel('design variable')
    title('gradients')
    grad
    view(0,90)
    subplot(2,2,2)
    grads=vertcat(supportOptim.hist(:).prevDir);
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
    grads=vertcat(supportOptim.hist(:).gradfk);
    gradsm1=vertcat(supportOptim.hist(:).gradfkm1);
    dot(grads,gradsm1,2);
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
    for ii=1:numel(postDir)
        
        dirDat(ii).path=postDir{ii};
        [dirDat(ii).case,dirDat(ii).supportOptim,dirDat(ii).diffStep,dirDat(ii).paramoptim]=...
            ExtractParamOptim(postDir{ii});
        
        [h]=PlotGradients(dirDat(ii).paramoptim);
        hgsave(h,[postDir{ii},filesep,'GradHistory_',h.Name])
        close(h)
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
    if isempty(regexp(MatName{1},'partial','once'))
        load(MatPath{1})
    else
        load(MatPath{1})
        
    end
    
    if exist('paramoptim','var')
        supportOptim=paramoptim.optim.supportOptim;
        caseName=paramoptim.general.optimCase;
        diffStep=paroptim.optim.CG.minDiffStep;
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