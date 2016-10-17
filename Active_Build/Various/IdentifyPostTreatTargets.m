function [T,postDir]=IdentifyPostTreatTargets(rootDir,maxRecurs,latestDate)
    
    [childFolder]=ExploreFolderTree(rootDir,maxRecurs);
    %[childFolder]=CheckDateEdit(childFolder,latestDate);
    [resDir]=FindResDirFolder(childFolder);
    [postDir]=NeedsPostTreat(resDir);
    
    [T]=PostTreatFiles(postDir);
    
end

function [oldDir]=CheckDateEdit(newDir,latestDate)
    
    latestDate=datenum(latestDate);
    isOldDir=false(size(newDir));
    for ii=numel(newDir)
        
        c=dir(newDir{ii});
        isOldDir=datenum(c(1).date)<latestDate;
    end
    
    oldDir=newDir(isOldDir);
end

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
            isempty(FindDir(resDir{ii},'Optimisation_',false));
        
    end
    postDir=resDir(isPost);
end

function [T,errDir]=PostTreatFiles(postDir)
    isError=false(size(postDir));
    fprintf('\n%i folders to be postreated \n',numel(postDir))
    for ii=1:numel(postDir)
        fprintf('%s\n',postDir{ii})
    end
    fprintf('\n')
    parfor ii=1:numel(postDir)
        
        try
            [optimstruct]=ExtractOptimStruct(postDir{ii});
            
            [T{ii}]=CaptureOutPostTreatIncomplete(postDir{ii},optimstruct)
            fprintf('   Folder %i - done ',ii)
        catch MEIdentif
            T{ii}=MEIdentif;
            isError(ii)=true;
            fprintf('   Folder %i - error ',ii)
        end
        
    end
    fprintf('\n')
    errDir=postDir(isError);
    if numel(errDir)>0
        warning('Errors detected, %i folders affected. See list below',sum(isError));
        for ii=1:numel(errDir)
            fprintf('%s\n',errDir{ii})
        end
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