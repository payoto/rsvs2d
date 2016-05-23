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
    
    varExtract={'direction','optimMethod'};
    [direction,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    if numel(nIter)>1
        nIterStart=nIter(1);
        nIterEnd=nIter(2);
        nIter=nIter(2)-nIter(1);
    end
    
    iterstruct=struct([]);
    iterstruct=repmat(iterstruct,[nIter,1]);
    
    for ii=nIterStart:nIterEnd
        
        [iterPath,~]=FindDir(pathStr,['iteration_',int2str(ii)],true);
        
        datName=['population_iteration_',int2str(ii)];
        [returnPath,~]=FindDir(iterPath{1},datName,false);
        load(returnPath{1},'population');
        iterstruct(ii-nIterStart+1).population=population;
        
    end
    
    if TestGreed(optimMethod)
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
    
end

function [isgreedy]=TestGreed(optimMethod)
    
   isgreedy=true;
   switch optimMethod
       case 'DE'
       case 'DEtan'
       case 'conjgrad'
           isgreedy=false;
       case 'conjgradls'
           isgreedy=false;
       otherwise
           error('Need to say if optimisation method is greedy or not')
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
