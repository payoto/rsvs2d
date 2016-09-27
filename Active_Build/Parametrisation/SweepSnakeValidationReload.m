
function [procdat,parampreset,sweepopts]=SweepSnakeValidationReload(sweepCaseStr,rootDir)
    
    include_Utilities
    include_Validation
    
    sweepopts=eval(sweepCaseStr{1});
    param=eval(sweepCaseStr{2});
    
    [sweepopts]=GenerateSweepArrays(sweepopts);
    [parampreset]=GenerateSweptParam(sweepopts,param);
    procdat=repmat(ProcDatTemplate,[length(parampreset),12]);
    
    for ii=1:length(parampreset)
        
        [procdatTemp]=procdatfindProcDat(rootDir,[[sweepCaseStr{:}],'_',int2str(ii)]);
        nDat(ii)=length(procdatTemp);
        procdat(ii,1:nDat(ii))=procdatTemp;
    end
    procdat=procdat(:,1:max(nDat));
    OutputSweepData([sweepCaseStr{:}],procdat,sweepopts,parampreset)
    
end

function [procdat]=procdatfindProcDat(rootDir,str)
    [returnPath,returnName]=FindDir(rootDir,str,true);
    [returnPath,returnName]=FindDir(returnPath{1},'Procdat',false);
    load(returnPath{1})
    
end

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
        disp('FindDir Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end


%% Output

function []=OutputSweepData(sweepCaseStr,procdat,sweepopts,parampreset)
    
    param=structInputVar('CurrentValidation');
    varExtract={'resultRoot','archiveName'};
    [resultRoot,archiveName]=ExtractVariables(varExtract,param);
    
    [marker,t]=GenerateResultMarker(sweepCaseStr);
    [writeDirectory]=GenerateValSweepDirectory(marker,resultRoot,archiveName,t);
    
    fileName=[writeDirectory,filesep,'Procdat_',marker,'.mat'];
    save(fileName,'procdat','sweepopts','parampreset');
    
    paramCell{1}=['# ',sweepCaseStr,' - ',datestr(t),' - ',marker];
    paramCell{2}=['# '];
    
    [tabLines]=GenerateSweepTable(procdat,sweepopts,parampreset);
    
    OutputSweepTableResult([paramCell,tabLines],writeDirectory,marker)
    
end

function [resultDirectory]=GenerateValSweepDirectory(marker,resultRoot,...
        archiveName,t)
    if ~exist('t','var'),t=now;end
    dateSubFolders=['Archive_',datestr(now,'yyyy_mm'),'\Day_',datestr(t,29)];
    resultDirectory=[resultRoot,filesep,archiveName,filesep,dateSubFolders,...
        filesep,'SWEEP_',marker];
    
    resultDirectory=MakePathCompliant(resultDirectory);
    
    mkdir(resultDirectory)
end


function []=OutputSweepTableResult(tableCell,writeDirectory,marker)
    
    fileName=['sweeptable_',marker,'.csv'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
    WriteToFile(tableCell,FID)
    
    fclose(FID);
end

function [tabLines]=GenerateSweepTable(procdat,sweepopts,parampreset)
    
    varNames={sweepopts(:).varname};
    
    tabLines{size(procdat,1)+1}='';
    
    for ii=1:size(procdat,1)
        
        [sumLines]=GenerateSummaryLine(procdat(ii,:));
        [paramSubstring]=BuildParamSubstring(varNames,parampreset(ii));
        if ii==1
            tabLines{1}=[paramSubstring{1},',',sumLines{1}];
        end
        tabLines{ii+1}=[paramSubstring{2},',',sumLines{2}];
    end
    
end

function [paramSubstring]=BuildParamSubstring(varNames,parampreset)
    
    paramSubstring{1}='';
    paramSubstring{2}='';
    separator='';
    for ii=1:length(varNames)
        paramVal=ExtractVariables(varNames(ii),parampreset);
        paramSubstring{1}=[paramSubstring{1},separator,varNames{ii}];
        paramSubstring{2}=[paramSubstring{2},separator,ProcesstoString(paramVal)];
        separator=' , ';
    end
    
    
end
%% Generate Swept Param

function [sweepopts]=GenerateSweepArrays(sweepopts)
    
    for ii=1:length(sweepopts)
        if isempty(sweepopts(ii).array)
            switch sweepopts(ii).steps{2}
                case 'lin'
                    sweepopts(ii).array=linspace(sweepopts(ii).range(1),...
                        sweepopts(ii).range(2),sweepopts(ii).steps{1});
                case 'log'
                    sweepopts(ii).array=logspace(sweepopts(ii).range(1),...
                        sweepopts(ii).range(2),sweepopts(ii).steps{1});
            end
        end
    end
    
    
end

function [param]=GenerateSweptParam(sweepopts,param)
    
    nSteps=length([sweepopts(:).array]);
    
    
    param.structdat=GetStructureData(param);
    param=repmat(param,[nSteps,1]);
    
    kk=1;
    
    for ii=1:length(sweepopts)
        for jj=1:length(sweepopts(ii).array)
            if ~iscell(sweepopts(ii).array)
            param(kk)=SetVariables({sweepopts(ii).varname},...
                {sweepopts(ii).array(jj)},param(kk));
            else
                param(kk)=SetVariables({sweepopts(ii).varname},...
                {sweepopts(ii).array{jj}},param(kk));
            end
            kk=kk+1;
        end
    end
    
    
end

%% Sweep Cases
%% Sweep Cases

function [sweepopts]=lEps_arrTol()
    
    kk=1;
    sweepopts(kk).varname='lengthEpsilon';
    sweepopts(kk).array=[];
    sweepopts(kk).range=[-4 -8];
    sweepopts(kk).steps={[5],'log'};
    kk=kk+1;
    
    sweepopts(kk).varname='arrivalTolerance';
    sweepopts(kk).array=[];
    sweepopts(kk).range=[-4 -1];
    sweepopts(kk).steps={[4],'log'};
    kk=kk+1;
        
    
end

function [sweepopts]=lEps_InitPos()
    
    kk=1;
    sweepopts(kk).varname='lengthEpsilon';
    sweepopts(kk).array=[];
    sweepopts(kk).range=[-4 -6];
    sweepopts(kk).steps={[5],'log'};
    kk=kk+1;
    
    sweepopts(kk).varname='snaxInitPos';
    sweepopts(kk).array=[];
    sweepopts(kk).range=[-4 -6];
    sweepopts(kk).steps={[5],'log'};
    kk=kk+1;
        
    
end

function [sweepopts]=test1()
    
    kk=1;
    sweepopts(kk).varname='lengthEpsilon';
    sweepopts(kk).array=[1e-4 1e-8];
    sweepopts(kk).range=[-4 -8];
    sweepopts(kk).steps={[2],'log'};
    kk=kk+1;
    
        
    
end

%%
function [param]=ActiveParameters()
    
    
    % Note Files
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','Opimisation','VALIDATION','SQP','Profile Length'};
    param.results.archiveName='ParamValidation';
    
    % Local optimum avoidance params
    param.snakes.step.mergeTopo=true;

    param.snakes.force.typeSmear='length';
    param.snakes.step.arrivalTolerance=1e-1;
    param.snakes.force.lengthEpsilon=1e-5;
    param.snakes.step.snaxInitPos=1e-5;
    param.snakes.step.convCheckRate=20;
    param.snakes.step.convCheckRange=15;
    param.snakes.step.convDistance=200;
    param.snakes.step.fillLooseStep=5;
    param.snakes.step.fillLooseCut=1e-3;
    
    % Default stepping params for validation (some cases might need more)
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    % Need to sort out the domain sizes to be always square?
end

function [param]=ActiveParameters2()
    
    
    % Note Files
    param.results.noteFiles={'CurrentBuild'};
    param.results.tags={'snakes','Opimisation','VALIDATION','SQP','Profile Length'};
    param.results.archiveName='ParamValidation';
    
    % Local optimum avoidance params
    param.snakes.step.mergeTopo=true;

    param.snakes.force.typeSmear='d';
    param.snakes.step.arrivalTolerance=1e-1;
    param.snakes.force.lengthEpsilon=1e-5;
    param.snakes.step.snaxInitPos=1e-5;
    param.snakes.step.convCheckRate=20;
    param.snakes.step.convCheckRange=15;
    param.snakes.step.convDistance=200;
    param.snakes.step.fillLooseStep=5;
    param.snakes.step.fillLooseCut=1e-3;
    
    % Default stepping params for validation (some cases might need more)
    param.snakes.step.snakesSteps=150;
    param.snakes.refine.refineGrid=4;
    param.snakes.refine.typeRefine='all';
    
    % Need to sort out the domain sizes to be always square?
end



