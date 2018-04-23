function [optimstruct]=RecoverASOFlow(pathStr,iterN)
    
    optPath=FindDir(pathStr,'OptimRes',0);
    paramPath=FindDir(pathStr,'FinalParam',0);
    optimstruct=load(optPath{1});
    optimstruct=optimstruct.optimstruct;
    paramoptim=load(paramPath{1});
    paramoptim=paramoptim.paramoptim;
    
    for ii=iterN
        for jj=1:numel(optimstruct(ii).population)
            try
                optimstruct(ii).population(jj)=RecoverASORun(optimstruct(ii).population(jj),paramoptim);
                optimstruct(ii).population(jj).constraint=1;
            catch MEid
                disp([int2str(ii),':',int2str(jj),' failed'])
                MEid.getReport
            end
        end
%         tecFilePath=FindDir(pathStr,'Tec360PLT_',0);
%         tecFilePath=MakePathCompliant(tecFilePath{1});
%         if ExtractVariables({'useSnake'},paramoptim)
%             axisRatio=ExtractVariables({'axisRatio'},paramoptim.parametrisation);
%             ConcatenateTecplotFile(pathStr,tecFilePath,axisRatio)
%         end
    end
    
    outinfo=ReconstructOutinfo(optimstruct);
    nIter=numel(optimstruct);
    while isempty(optimstruct(nIter).population(1).objective)
        nIter=nIter-1;
    end
        
    OptimisationOutput('final',paramoptim,outinfo,optimstruct(1:nIter));
end


function member=RecoverASORun(member,paramoptim)
    
    ASOresult=ASOPart([member.location,filesep,'workspace.mat']);
    loopPath=FindDir(member.location,'restart_',0);
    
    load(loopPath{1},'loop');
    
    % ---------------------
    % Organise Outputs
    
    %[~,areaAdd]=LengthArea(paramoptim,member,loop);
    obj=ASOresult.flow;
    objValue=obj.CD;
    additional=obj;
    
    
    
    % ---------------------
    % Calculate Fill Movement
    
    varExtract={'axisRatio'};
    [axisRatio]=ExtractVariables(varExtract,paramoptim.parametrisation);
    varExtract={'typeLoop'};
    [typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    varExtract={'asoCase','asoPath','nMach','asoReturnFillChange',...
        'desVarConstr','desVarVal','su2ProcSec','asoProcSec','snoptIter'};
    [asoCase,asoPath,nMach,asoReturnFillChange,desVarConstr,desVarVal,...
        su2ProcSec,asoProcSec,snoptIter]=ExtractVariables(varExtract,paramoptim);
    %error('Returning Fill delta not coded yet')
    if isempty(ASOresult.loop)
        ASOresult.loop=struct('coord',[]);
        for ii=1:numel(loop)
            ASOresult.loop(ii).coord=(loop(ii).(typeLoop));
        end
        warning('ASOresult is not returning the loop')
    elseif ~isfield(ASOresult.loop,'ASOResult') && isstruct(ASOresult.loop)
        [ASOresult.loop.coord]=deal(ASOresult.loop.(typeLoop));
    else
        [ASOresult.loop.coord]=deal(ASOresult.loop.ASOResult);
        
    end
    for ii=1:numel(ASOresult.loop)
        ASOresult.loop(ii).coord(:,2)=ASOresult.loop(ii).coord(:,2)/axisRatio;
    end
    [~,areaAdd]=LengthArea(paramoptim,member,{ASOresult.loop.coord});
    if asoReturnFillChange
        
        [fill,~]=LoopToFill(ASOresult.loop,baseGrid);
        additional.filldelta=fill-member.fill;
    else
        additional.filldelta=member.fill-member.fill;
    end
    
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
    
    
    fieldsAdd=fieldnames(additional);
    filldeltaInd=find(~cellfun(@isempty,regexp(fieldsAdd,'filldelta')));
    if ~isempty(filldeltaInd)
        member.fill=member.fill+additional.filldelta;
    end
    % Additional Data to be saved
    fieldsAdd(filldeltaInd)=[];
    for jj=1:numel(fieldsAdd)
        member.additional.(fieldsAdd{jj})=additional.(fieldsAdd{jj});
    end
    member.objective=objValue;
end


function [result]=ASOPart(pathStr)
    load(pathStr)
    if ~flagNoASO
        xstar = optimiser.results.x;
        geom = surface.deParam(xstar);
        loops = geom.getLoops();
        for i=1:length(loops)
            loop(i).ASOResult = geom.vertices(loops{i}(1:end-1,:),:);
        end %for
    else
        loop=[];
    end
    vars = {'CL','CD','CM','directResidual','directNIter','adjointResidual',...
        'adjointNIter','majorIt','minorIt'};
    for i=1:length(vars)-2
        result.flow.(vars{i}) = solver.(vars{i});
    end %for
    for i=length(vars)-1:length(vars)
        result.flow.(vars{i}) = log.(vars{i});
    end %for
    result.loop = loop;
    result.log = log;
    
end


function [objValue,additional]=LengthArea(paramoptim,member,loop)
    if isstruct(loop)
        for ii=1:length(loop)
            
            [xMin(ii),xMax(ii),t(ii),L(ii),A(ii)]=...
                ClosedLoopProperties(loop(ii).snaxel.coord(1:end-1,:));
            
        end
    elseif iscell(loop)
        for ii=1:length(loop)
            
            [xMin(ii),xMax(ii),t(ii),L(ii),A(ii)]=...
                ClosedLoopProperties(loop{ii});
            
        end
    end
    objValue=sum(A)/sum(L);
    
    additional.A=sum(A);
    additional.L=sum(L);
    additional.t=sum(t);
    additional.c=max(xMax)-min(xMin);
    additional.tc=additional.t/additional.c;
end

function []=ConcatenateTecplotFile(iterDir,tecFilePath,axisRatio)
    [profPaths]=FindProfile(iterDir,axisRatio);
    
    compType=computer;
    
    if strcmp(compType(1:2),'PC')
        for ii=1:length(profPaths)
            system(['type "',profPaths{ii},'" >> "',tecFilePath,'"']);
        end
    else
        for ii=1:length(profPaths)
            system(['cat ''',profPaths{ii},''' >> ''',tecFilePath,'''']);
        end
    end
    
    
end

function [profPaths]=FindProfile(iterDir,axisRatio)
    
    [returnPath]=FindDir(iterDir,'profile',true);
    
    for ii=1:length(returnPath)
        
        profPaths(ii)=FindDir(returnPath{ii},'tecsubfile',false);
        try
            surfPlt=[returnPath{ii},filesep,'run',filesep,'surface.plt'];
            if exist(surfPlt,'file')
                MergeTecSubfile(profPaths{ii},surfPlt,axisRatio);

            end
        catch ME
            ME.getReport
        end
    end
    
end