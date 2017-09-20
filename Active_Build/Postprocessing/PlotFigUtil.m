function [pltPaths,h]=PlotFigUtil(rootDir,dir,nameStr)
    
    
    [childFolder]=ExploreFolderTree(rootDir);
    
%    FindDir([minIterPos,filesep,'CFD'],'Tec360plt_',false)
    [pltPaths]=FindMATFiles(childFolder);
    h=PreplotFiles(pltPaths,dir,nameStr);
end

function [h]=PreplotFiles(pltPaths,dir,nameStr)
   
    for ii=1:length(pltPaths)
        
        [h(ii)]=PlotFiles(pltPaths{ii},dir,nameStr);
        
    end
    
end

function [h]=PlotFiles(matPath,dir,nameStr)
    
    if ~exist('nameStr','var'),
        nameStr='Optimisation History';
    else
        nameStr=eval(nameStr);
        nameStr=nameStr{end};
    end
    optimstruct=[];
    load(matPath,'optimstruct');
    ii=1;
    while (ii<=numel(optimstruct) && ~isempty(optimstruct(min(ii,numel(optimstruct))).population(1).objective))
        ii=ii+1;
    end
    h=0;
    if ii>3
        [h]=OptimHistory(true,optimstruct(1:ii-1),0,1000,dir);
        nameStr=regexprep(nameStr,'\.mat','');
        h.Name=nameStr;
        figPath=regexprep(matPath,'OptimRes_','Optimisation_');
        figPath=regexprep(figPath,'.mat','.fig');
        hgsave(h,figPath);
    end 
    
    
    
end

function [pltPaths]=FindMATFiles(childFolder)
    kk=1;
    for ii=1:length(childFolder)
        intermPath=FindDir(childFolder{ii},'OptimRes_',false);
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