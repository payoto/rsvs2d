function [objValue,additional]=ASOShape(paramoptim,member,loop,baseGrid)
    % Function for the interface of the RSVS shape optimisation framework
    % and the B-Spline and Smoothness constraint ASO
    %
    
    
    
    % ---------------------
    % Generate convHull boundary
    
    % ---------------------
    % Interface parameters
    thisworker = getCurrentWorker;
    varExtract={'workerList','machineList'};
    [workerList,machineList]=ExtractVariables(varExtract,paramoptim);
    varExtract={'axisRatio','typeLoop'};
    [axisRatio,typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    if ~isempty(thisworker)
        currentMachineFile=machineList(thisworker.ProcessId==workerList);
    else
        try
        currentMachineFile=machineList(1);
        catch
            machineList=WriteMachineFiles(1,'test');
            currentMachineFile=machineList(1);
        end
    end
    optimDirectory=member.location;
    
    
    % ---------------------
    % Match ASO and Current Optimisation
    varExtract={'asoCase','asoPath','nMach','asoReturnFillChange',...
        'desVarConstr','desVarVal','su2ProcSec','asoProcSec','snoptIter',...
        'aeroName','profileComp'};
    [asoCase,asoPath,nMach,asoReturnFillChange,desVarConstr,desVarVal,...
        su2ProcSec,asoProcSec,snoptIter,aeroName,profileComp]=ExtractVariables(varExtract,paramoptim);
    
    addpath(MakePathCompliant(asoPath));
    addpath(MakePathCompliant([asoPath,filesep,'matlab-snopt']));
    
    
    ASOOptions = asoCase();
    
    ASOOptions.snopt.wcLimit = ceil(asoProcSec/currentMachineFile.slots);
    ASOOptions.snopt.maxIter = snoptIter;
    
   
    targetLoop=NacaMultiTopo(1000,aeroName);
    ASOOptions.shapeMatching.initialLoops=CloseCellLoops({loop.(typeLoop)});
    ASOOptions.shapeMatching.targetLoops=CloseCellLoops({targetLoop.coord});
    
    
    
    ASOOptions.shapeMatching.Fcn=@(targetLoops,initialLoops) ...
        InverseDesign_ErrorTopoAD(profileComp,initialLoops,targetLoops);
    
    % ASOOptions.solver.mpiOpts=['--hostfile "','mpihostfile','"  --oversubscribe'];
    
    % Other Options
    ASOOptions.structdat=GetStructureData(ASOOptions);
    paramoveride=paramoptim.obj.aso.paramoveride;
    fieldsOverride=fieldnames(paramoveride);
    
    for ii=1:numel(fieldsOverride)
        ASOOptions=SetVariables(fieldsOverride(ii),...
            {paramoveride.(fieldsOverride{ii})},ASOOptions);
    end
    
    % Call the correct constraints
    for ii=1:numel(desVarConstr)
        switch desVarConstr{ii}
            case 'subdiv_ebasis_cross'
                ASOOptions.problemargin=[ASOOptions.problemargin,...
                    {'subdiv_ebasis_cross',{}}];
                
            case 'subdiv_ebasis_lim'
                ASOOptions.problemargin=[ASOOptions.problemargin,...
                    {'subdiv_ebasis',{}}];
            case 'ValVolFrac'
            case 'MinSumVolFrac'
            case 'MaxSumVolFrac'
            case 'MinValVolFrac'
                ASOOptions.problemargin=[ASOOptions.problemargin,...
                    {'area_gt',desVarVal(ii)}];
            case 'LocalVolFrac_min'
            case 'LocalVolFrac_equal'
            case 'LocalVolFrac_max'
        end
    end
    
    ASOOptions=rmfield(ASOOptions,'structdat');
    % --------------------------
    % Call ASO
    
    origDir=cd(optimDirectory);
    optimDirectory=['.'];
    ASOresult = ShapeMatch(optimDirectory,ASOOptions);
    optimDirectory=cd(origDir);
    
    
    % ---------------------
    % Organise Outputs
    
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    obj=ASOresult.flow;
    objValue=obj.CD;
    additional=obj;
   
    
    
    % ---------------------
    % Calculate Fill Movement
    
    varExtract={'typeLoop'};
    [typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    
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
end

function [cellLoop]=CloseCellLoops(cellLoop)
    
    for ii=1:numel(cellLoop)
        if any(cellLoop{ii}(end,:)~=cellLoop{ii}(1,:))
            cellLoop{ii}(end+1,:)=cellLoop{ii}(1,:);
        end
    end
    
end

function [loopconstr]=PushConvHull(loop,typeLoop)
    loopconstr=repmat(struct('coord',zeros([0,2]),'coord',zeros([0,2]),...
        'coord',zeros([0,2])),size(loop));
    paramspline.splineCase='convhulltri';
    for ii=1:numel(loop)
        loopconstr(ii).coord=loop(ii).(typeLoop);
        loopconstr(ii).constraint=loopconstr(ii).coord(...
            convhull(loopconstr(ii).coord(:,1),loopconstr(ii).coord(:,2)),:);
        %         paramspline.splineCase=min(max(
        loopconstr(ii).constraint=ResampleSpline
        
    end
    
    
end

function []=MeshTriangleSU2()
    
    boundaryLoc=member.location;
    SU2Flow_Handler(paramoptim,boundaryLoc);
    copyfile([boundaryLoc,filesep,'SU2CFD',filesep,'triangularmesh.su2'],...
        [boundaryLoc,filesep,'mesh.su2'])
    rmdir([boundaryLoc,filesep,'SU2CFD'],'s')
    
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