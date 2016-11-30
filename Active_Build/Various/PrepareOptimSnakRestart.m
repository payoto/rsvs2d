function [fOut]=PrepareOptimSnakRestart(pathStr,numIter,numProf)
    
    
    fileName=[pathStr,filesep,'iteration_',int2str(numIter),filesep,'profile_',...
        int2str(numProf),filesep,'restart_',int2str(numIter),'_',int2str(numProf)];
    
    load(fileName)
    snakrestart=restartsnak;
    
    fileName=[pathStr,filesep,'iteration_',int2str(0),filesep,'profile_',...
        int2str(0),filesep,'restart_',int2str(0),'_',int2str(0)];
    load(fileName)
    gridrefined=grid.refined;
    connectstructinfo=grid.connec;
    unstructReshape=grid.base;
    snakrestart.snaxel=restartsnak.snaxel;
    snakrestart.insideContourInfo=restartsnak.insideContourInfo;
    [returnPath,~]=FindDir(pathStr,'FinalParam',0);
    load(returnPath{1})
    param=paramoptim.parametrisation;
    param=SetVariables({'snakData'},{'all'},param);
    
    fOut=['OptimSnakRestart_',int2str(numIter),'_',int2str(numProf)];
    save(fOut,'gridrefined','snakrestart','param','connectstructinfo','unstructReshape');
    
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