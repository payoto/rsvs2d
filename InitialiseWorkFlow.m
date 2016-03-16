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


function []=InitialiseWorkFlow()
    
    singleFolder={'Active_Build'};
    rootTreeFolders={'Active_Build\Input_Variables','Active_Build\Include_Files','Active_Build\Various',...
        'Active_Build\Velocity','MEX_Function_Directory\MEX_Executables','Automated_Function_Directory'};
    
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
        addFolders{ii}=[cd,'\',folders{ii}];
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
                branchDir{jj}=[rootDir{ii},'\',dirNames{jj}];
            end

            [addSubFolders]=ExploreFolderTree(branchDir);
            addFolders=[addFolders,addSubFolders];
        end
    end
    
end

function []=AddFoldersToPath(addFolders)
    % adds a set of paths to the active path
    
    for ii=1:length(addFolders)
        addpath(addFolders{ii})
    end
    
    
end


