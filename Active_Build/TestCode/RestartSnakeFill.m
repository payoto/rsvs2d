% Function which restarts a failed snaking process from the appropriate
% data.

function [grid,loop,restartsnak,snakSave,newFill]=RestartSnakeFill(...
        optimstruct,nIter,nProf,rootFolder,restartType)
    % Runs if 3rd input of ExecuteOptimisation exists.
    [grid,loop,restartsnak,snakSave]=ExtractRestartPaths(optimstruct,nIter,restartType,rootFolder);
    % get grid information
    
    % run as in ExecuteOptimisation
    if isstruct(optimstruct)
        newFill=optimstruct(nIter).population(nProf).fill;
    else
        newFill=optimstruct{2};
    end
    
end

function [grid,loop,restartsnak,snakSave]=ExtractRestartPaths(optimstruct,nIter,restartType,rootFolder)
    
    
    pathGrid=[rootFolder,filesep,'iteration_0',filesep,'profile_0',filesep];
    if strcmp(restartType,'iter0')
        pathRestart=pathGrid;
    elseif strcmp(restartType,'sensitivity')
        if isstruct(optimstruct)
            pathRestart=optimstruct(nIter).population(1).location;
        else
            pathRestart=optimstruct{1};
        end
    else
        error('Invalid type')
    end
    
    [returnPath,returnName]=FindDir(pathGrid,'restart',0);
    load(returnPath{1})
    [returnPath,returnName]=FindDir(pathRestart,'restart',0);
    load(returnPath{1})
   
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
        disp('FindDir Could not find requested item')
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


function [newGrid,newRefGrid]=ReFracGrids(baseGrid,refinedGrid,...
        connectstructinfo,newFracs)
    
    activeCell=true(size(baseGrid.cell));
    activeInd=[baseGrid.cell((activeCell)).index];
    
    connecInd=[connectstructinfo.cell(:).old];
    activConnecSub=FindObjNum([],activeInd,connecInd);
    activeCellSub=find(activeCell);
    refCellInd=[refinedGrid.cell(:).index];
    
    if numel(newFracs)~=numel(activeCellSub)
        error('Fill and Active Set do not match in size')
    end
    if sum(abs(newFracs))==0
        newFracs(round(numel(newFracs)/2))=1e-3;
    end
    
    newGrid=baseGrid;
    newRefGrid=refinedGrid;
    
    for ii=1:length(activeCellSub)
        newGrid.cell(activeCellSub(ii)).fill=newFracs(ii);
        newSub=FindObjNum([],[connectstructinfo.cell(activConnecSub(ii)).new],refCellInd);
        [newRefGrid.cell(newSub).fill]=deal(newFracs(ii));
    end
    
end

