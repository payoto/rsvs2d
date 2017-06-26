% Function for the postreatment of incomplete optimisation run
% The goal is to reconstruct enough information to take the run through 

function [outinfo,paramoptim,iterstruct]=PostTreatIncomplete(pathStr,nIter,iterstruct)
    include_Utilities
    include_PostProcessing
    include_Optimisation
    
    % Reconstruct outinfo
%     outinfo.rootDir=pathStr;
%     [outinfo.tOutput,outinfo.marker]=FindTime(pathStr);
    [outinfo]=ReconstructOutinfo(iterstruct);
    % iterstruct paramoptim
    [paramoptim2]=ReconstructParameter(pathStr,outinfo(1).marker);
    try 
        paramoptim=StructOptimParam(ExtractVariables({'optimCase'},paramoptim2));
        outVars={paramoptim2.structdat(:).vars.name};
        for ii=numel(outVars):-1:1
            [valVars{ii}]=ExtractVariables(outVars(ii),paramoptim2);
        end
        for ii=1:numel(outVars)
            try
            paramoptim=SetVariables(outVars(ii),valVars(ii),paramoptim);
            catch
            end
        end
    catch MEId
        disp(MEId.getReport)
        [paramoptim]=ReconstructParameter(pathStr,outinfo(1).marker);
    end
    % Reconstruct iterstruct
    if nargin<3
        [iterstruct]=ReconstructIterationStructure(pathStr,nIter,paramoptim);
    end
    
    
    OptimisationOutput('final',paramoptim,outinfo,iterstruct);
end


function [outinfo]=ReconstructOutinfo(optimstruct)
    
    allRootDir=repmat(struct('rootDir',''),[1 numel(optimstruct)]);
    rmR=[];
    for ii=1:numel(optimstruct)
        try
            kk=1;
            while isempty(optimstruct(ii).population(kk).location)
                kk=kk+1;
            end
            allRootDir(ii).rootDir=regexprep(optimstruct(ii).population(kk).location,'iteration.*$','');
        catch
            
        
            rmR=[ii];
        end
    end
    allRootDir(rmR)=[];
    [allRootDir,uniqueRootDir]=IdentifyUniqueOptions(allRootDir);
    outinfo=repmat(struct('rootDir','','tOutput',[],'marker',''),[1 numel(uniqueRootDir{1})]);
    [outinfo(:).rootDir]=deal(uniqueRootDir{1}{:});
    for ii=1:numel(outinfo)
        [outinfo(ii).rootDir,outinfo(ii).tOutput,outinfo(ii).marker]=FindTime2(outinfo(ii).rootDir);
    end
    [~,sortOrd]=sort([outinfo(:).tOutput]);
    outinfo=outinfo(sortOrd);
end

function [refinestruct,patternList]=IdentifyUniqueOptions(refinestruct,jobfields)
    
    if nargin<2
        jobfields=fieldnames(refinestruct);
    end
    
    %
    
    n=numel(refinestruct);
    
    patternList=cell(0);
    for ii=1:numel(jobfields);
        patternList{ii}={};
        
        for jj=1:n
            optfound=false;
            
            for kk=1:numel(patternList{ii})
                test=false;
                if ischar(refinestruct(jj).(jobfields{ii}))
                    test=strcmp(patternList{ii}{kk},refinestruct(jj).(jobfields{ii}));
                elseif isnumeric(refinestruct(jj).(jobfields{ii})) || islogical(refinestruct(jj).(jobfields{ii}))
                    test=(patternList{ii}{kk}==refinestruct(jj).(jobfields{ii}));
                end
                if test
                    refinestruct(jj).([jobfields{ii},'num'])=kk;
                    optfound=true;
                    break
                end
            end
            if ~optfound
                patternList{ii}{end+1}=refinestruct(jj).(jobfields{ii});
                refinestruct(jj).([jobfields{ii},'num'])=numel(patternList{ii});
            end
        end
    end
    
    
    
    
end

function [t,marker]=FindTime(pathStr)
    
    [returnPath,returnName]=FindDir(pathStr,'Comments',false);
    
    fid=fopen(returnPath{1},'r');
    
    t=datenum(fgetl(fid));
    
    returnName=regexprep(returnName{1},'Comments_','');
    marker=regexprep(returnName,'.txt','');
    
end


function [pathStr,t,marker]=FindTime2(pathStr)
    
    splitPath=regexp(pathStr,'[\\,/]','split');
    pathStr=regexprep(pathStr,'[\\,/]$','');
    
    dirLoc=flip(find(~cellfun(@isempty,regexp(splitPath,'Dir'))));
    dirName=splitPath{dirLoc};
    splitDir=regexp(dirName,'_');
    t=datenum(dirName(splitDir(1)+1:splitDir(2)-1),'yyyy-dd-mmTHHMMSS');
    marker=dirName(splitDir(2)+1:end);
end

function [paramoptim]=ReconstructParameter(pathStr,marker)
    try
        [returnPath,returnName]=FindDir(pathStr,'FinalParam',0);
        load(returnPath{1})
        if ~exist('paramoptim','var')
            error('Failure to load')
        end
    catch
        
        [rootParam,nameParam]=FindDir(pathStr,'param',0);
        isOptParam=cellfun(@isempty,regexp(nameParam,'_parametrisation'));
        matName=matlab.lang.makeValidName(regexprep(nameParam{isOptParam},'.dat',''));
        copyfile(rootParam{isOptParam},['.',filesep,matName,'.m']);
        eval(matName)
        paramoptim=param;
        paramoptim.structdat=GetStructureData(paramoptim);
        
        copyfile(rootParam{~isOptParam},['.',filesep,matName,'.m']);
        eval(matName)
        paramoptim.parametrisation=param;
        paramoptim.parametrisation.structdat=...
            GetStructureData(paramoptim.parametrisation);
    end
end

function [iterstruct]=ReconstructIterationStructure(pathStr,nIter,paramoptim)
    
    varExtract={'direction','optimMethod'};
    [direction,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    if numel(nIter)>1
        nIterStart=nIter(1);
        nIterEnd=nIter(2);
        nIter=nIter(2)-nIter(1)+1;
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
