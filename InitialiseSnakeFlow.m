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


<<<<<<< HEAD:InitialiseWorkFlow.m
function []=InitialiseWorkFlow()
    
    singleFolder={'Active_Build'};
    rootTreeFolders={'Active_Build\Input_Variables','Active_Build\Various','Active_Build\Velocity'};
=======
function []=InitialiseSnakeFlow()
    comStr=computer;
    if strcmp(comStr(1:2),'PC')
        
    else
        clear all
        setenv('TMP','/local/') 
    end
    
    singleFolder={''};
    rootTreeFolders={'Active_Build',['MEX_Function_Directory',filesep,'MEX_Executables']...
        ,'Automated_Function_Directory','JobScripts'};
>>>>>>> init_v2:InitialiseSnakeFlow.m
    
    [addSingleDir]=FormulateValidFolders(singleFolder);
    [rootTreeDirs]=FormulateValidFolders(rootTreeFolders);
    [branchesDir]=ExploreFolderTree(rootTreeDirs);
    
    addFolders=[addSingleDir,rootTreeDirs,branchesDir];
    
    AddFoldersToPath(addFolders)
end

function [addFolders]=FormulateValidFolders(folders)
    % adds a set of paths to the active path
    addFolders={};
    for ii=1:length(folders)
<<<<<<< HEAD:InitialiseWorkFlow.m
        addFolders{ii}=[cd,'\',folders{ii}];
=======
        addFolders{ii}=[cd,filesep,folders{ii}];
>>>>>>> init_v2:InitialiseSnakeFlow.m
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
<<<<<<< HEAD:InitialiseWorkFlow.m
                branchDir{jj}=[rootDir{ii},'\',dirNames{jj}];
=======
                branchDir{jj}=[rootDir{ii},filesep,dirNames{jj}];
>>>>>>> init_v2:InitialiseSnakeFlow.m
            end

            [addSubFolders]=ExploreFolderTree(branchDir);
            addFolders=[addFolders,addSubFolders];
        end
    end
    
end

function []=AddFoldersToPath(addFolders)
    % adds a set of paths to the active path
<<<<<<< HEAD:InitialiseWorkFlow.m
    
    for ii=1:length(addFolders)
        addpath(addFolders{ii})
    end
    
=======
    newPaths=addFolders{1};
    for ii=2:length(addFolders)
        newPaths=[newPaths,pathsep,addFolders{ii}];
    end
    addpath(newPaths);
>>>>>>> init_v2:InitialiseSnakeFlow.m
    
end


