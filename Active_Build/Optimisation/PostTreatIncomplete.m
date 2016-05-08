% Function for the postreatment of incomplete optimisation run
% The goal is to reconstruct enough information to take the run through 

function [outinfo,paramoptim,iterstruct]=PostTreatIncomplete(pathStr,nIter,iterstruct)
    include_Utilities
    include_PostProcessing
    
    % Reconstruct outinfo
    outinfo.rootDir=pathStr;
    [outinfo.tOutput,outinfo.marker]=FindTime(pathStr);
    % iterstruct paramoptim
    [paramoptim]=ReconstructParameter(pathStr,outinfo.marker);
    % Reconstruct iterstruct
    if nargin<3
        [iterstruct]=ReconstructIterationStructure(pathStr,nIter,paramoptim);
    end
    
    
    OptimisationOutput('final',paramoptim,outinfo,iterstruct);
end

function [t,marker]=FindTime(pathStr)
    
    [returnPath,returnName]=FindDir(pathStr,'Comments',false);
    
    fid=fopen(returnPath{1},'r');
    
    t=datenum(fgetl(fid));
    
    returnName=regexprep(returnName{1},'Comments_','');
    marker=regexprep(returnName,'.txt','');
    
end

function [paramoptim]=ReconstructParameter(pathStr,marker)
    
    fileName=['param_',marker];
    matName=matlab.lang.makeValidName(fileName);
    
    copyfile([pathStr,filesep,fileName,'.dat'],['.',filesep,matName,'.m']);
    eval(matName)
    paramoptim=param;
   paramoptim.structdat=GetStructureData(paramoptim);
   
   fileName=['param_',marker,'_parametrisation'];
    matName=matlab.lang.makeValidName(fileName);
    
    copyfile([pathStr,filesep,fileName,'.dat'],['.',filesep,matName,'.m']);
    eval(matName)
    paramoptim.parametrisation=param;
   paramoptim.parametrisation.structdat=...
       GetStructureData(paramoptim.parametrisation);
end

function [iterstruct]=ReconstructIterationStructure(pathStr,nIter,paramoptim)
    
    varExtract={'direction'};
    [direction]=ExtractVariables(varExtract,paramoptim);
    
    
    iterstruct=struct([]);
    iterstruct=repmat(iterstruct,[nIter,1]);
    
    for ii=1:nIter
        
        iterPath=[pathStr,filesep,'iteration_',int2str(ii)];
        datName=['population_iteration_',int2str(ii),'.mat'];
        
        load([iterPath,filesep,datName],'population');
        iterstruct(ii).population=population;
        
    end
    
    for ii=2:nIter
        iterPrecedent=[iterstruct(ii-1).population(:).objective];
        iterCurr=[iterstruct(ii).population(:).objective];
        switch direction
            case 'min'
                iterCopy=iterCurr>iterPrecedent;
            case 'max'
                iterCopy=iterCurr<iterPrecedent;
        end
        
        [iterstruct(ii).population(iterCopy)]=deal(iterstruct(ii-1).population(iterCopy));
       
        
    end
    
end

%% Standard Postreatment functions

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        warning('Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
      
    
    
end
