function [pltPaths]=PreplotUtil(rootDir)
    
    
    [childFolder]=ExploreFolderTree(rootDir);
    
%    FindDir([minIterPos,filesep,'CFD'],'Tec360plt_',false)
    [pltPaths]=FindPLTFiles(childFolder);
    PreplotFiles(pltPaths)
end

function []=PreplotFiles(pltPaths)
   
    for ii=1:length(pltPaths)
        
       launchcmd=['preplot "',pltPaths{ii},'" "',pltPaths{ii}(1:end-4),'_pre.plt"'];
       
       MaxSystemCommandInstances(launchcmd,4,'preplot')
        
    end
    
end

function [pltPaths]=FindPLTFiles(childFolder)
    kk=1;
    pltPaths=cell(0);
    for ii=1:length(childFolder)
        [intermPath,intermName]=FindDir(childFolder{ii},'Tec360plt_',false);
        isRem=false(size(intermName));
        for jj=1:numel(intermName)
            tempName=regexprep(intermName{jj},'\..*','');
            tempName=regexprep(tempName,'_pre$','');
            tempIsRem=~cellfun(@isempty,regexp(intermName,tempName));
            isRem([1:jj-1,jj+1:end])=isRem([1:jj-1,jj+1:end]) | tempIsRem([1:jj-1,jj+1:end]);
            isRem(jj)=isRem(jj) || any( tempIsRem([1:jj-1,jj+1:end]));
        end
        for jj=find(isRem)
            tempName=(regexprep(intermPath{jj},'\..*',''));
            if isempty(regexp(tempName,'_pre$','once'))
                c=dir([tempName,'.plt']);
                c(2)=dir([tempName,'_pre.plt']);
                isRem(jj)= ~(isempty(c(2).datenum) || (c(2).datenum<c(1).datenum));
            end
        end
        intermPath=intermPath(~isRem);
        if ~isempty(intermPath)
            pltPaths(kk:kk+numel(intermPath)-1)=intermPath;
            kk=kk+numel(intermPath);
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
