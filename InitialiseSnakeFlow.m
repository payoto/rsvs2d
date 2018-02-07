%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Initialise work Flow
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function []=InitialiseSnakeFlow(isRm)
    if nargin<1
        isRm=0;
    end
    
    comStr=computer;
    if strcmp(comStr(1:2),'PC')
        
    else
        clear all
        setenv('TMP','/local/')
        isRm=0;
    end
    
    singleFolder={''};
    rootTreeFolders={'Active_Build',['MEX_Function_Directory',filesep,'MEX_Executables']...
        ,'Automated_Function_Directory','JobScripts'};
    
    
    [addSingleDir]=FormulateValidFolders(singleFolder);
    [rootTreeDirs]=FormulateValidFolders(rootTreeFolders);
    [branchesDir]=ExploreFolderTree(rootTreeDirs);
    
    addFolders=[addSingleDir,rootTreeDirs,branchesDir];
    if isRm
        RmFoldersToPath(addFolders)
    else
        AddFoldersToPath(addFolders)
    end
end

function [addFolders]=FormulateValidFolders(folders)
    % adds a set of paths to the active path
    addFolders={};
    for ii=1:length(folders)
        
        addFolders{ii}=[cd,filesep,folders{ii}];
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

function []=AddFoldersToPath(addFolders)
    % adds a set of paths to the active path
    
    newPaths=addFolders{1};
    for ii=2:length(addFolders)
        newPaths=[newPaths,pathsep,addFolders{ii}];
    end
    addpath(newPaths);
end

function []=RmFoldersToPath(addFolders)
    % adds a set of paths to the active path
    
    newPaths=addFolders{1};
    for ii=2:length(addFolders)
        newPaths=[newPaths,pathsep,addFolders{ii}];
    end
    rmpath(newPaths);
end


