


function [procdat,parampreset,sweepopts]=SweepSnakeValidation(sweepCaseStr)
    
    include_Utilities
    include_Validation
    
    sweepopts=eval(sweepCaseStr);
    
    [sweepopts]=GenerateSweepArrays(sweepopts);
    [parampreset]=GenerateSweptParam(sweepopts);
    procdat=repmat(ProcDatTemplate,[length(parampreset),12]);
    
    for ii=1:length(parampreset)
        
        [procdatTemp,~]=SnakeValid([sweepCaseStr,'_',int2str(ii)],parampreset(ii));
        nDat(ii)=length(procdatTemp);
        procdat(ii,1:nDat(ii))=procdatTemp;
    end
    procdat=procdat(:,1:max(nDat));
    OutputSweepData(sweepCaseStr,procdat,sweepopts,parampreset)
    
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

function [param]=GenerateSweptParam(sweepopts)
    
    nSteps=length([sweepopts(:).array]);
    
    [param]=ActiveParameters();
    param.structdat=GetStructureData(param);
    param=repmat(param,[nSteps,1]);
    
    kk=1;
    
    for ii=1:length(sweepopts)
        for jj=1:length(sweepopts(ii).array)
            param(kk)=SetVariables({sweepopts(ii).varname},...
                {sweepopts(ii).array(jj)},param(kk));
            kk=kk+1;
        end
    end
    
    
end

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
    param.snakes.step.arrivalTolerance=10e-2;
    param.snakes.force.lengthEpsilon=1e-6;
    param.snakes.step.snaxInitPos=10*param.snakes.force.lengthEpsilon;
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






