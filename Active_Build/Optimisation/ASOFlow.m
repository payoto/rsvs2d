function [objValue,additional]=ASOFlow(paramoptim,member,loop,baseGrid)
    % Function for the interface of the RSVS shape optimisation framework
    % and the B-Spline and Smoothness constraint ASO
    %
    
    % ---------------------
    % Generate Mesh
    boundaryLoc=member.location;
    SU2Flow_Handler(paramoptim,boundaryLoc);
    copyfile([boundaryLoc,filesep,'SU2CFD',filesep,'triangularmesh.su2'],...
        [boundaryLoc,filesep,'mesh.su2'])
    rmdir([boundaryLoc,filesep,'SU2CFD'],'s')
    
    % ---------------------
    % Generate convHull boundary
    
    % ---------------------
    % Interface parameters
    thisworker = getCurrentWorker;
    varExtract={'workerList','machineList'};
    [workerList,machineList]=ExtractVariables(varExtract,paramoptim);
    if ~isempty(thisworker)
        currentMachineFile=machineList(thisworker.ProcessId==workerList);
    else
        currentMachineFile=machineList(1);
    end
    optimDirectory=boundaryLoc;
    
    
    % ---------------------
    % Call ASO
    varExtract={'asoCase','asoPath','nMach','asoReturnFillChange'};
    [asoCase,asoPath,nMach,asoReturnFillChange]=ExtractVariables(varExtract,paramoptim);
    addpath(MakePathCompliant(asoPath));
    addpath(MakePathCompliant([asoPath,filesep,'matlab-snopt']));
    copyfile([MakePathCompliant(asoPath),filesep,'templatesu2.cfg'],...
        [optimDirectory,filesep,'su2.cfg'])
    
    ASOOptions = asoCase();
    ASOOptions.solver.mach = nMach;
    ASOOptions.solver.np=currentMachineFile.slots;
    ASOOptions.problemargin={'area_gt'};
    
    copyfile(currentMachineFile.file,[optimDirectory,filesep,'mpihostfile'])
    ASOOptions.solver.mpiOpts=['--hostfile "','mpihostfile','"'];
    
    origDir=cd(optimDirectory);
    optimDirectory=['.'];
    ASOresult = ASO(optimDirectory,ASOOptions);
    % meshfile: SU2.SU2toPLT(path2Mesh)
    % meshfile: !SU2_SOL run/su2.cfg
    optimDirectory=cd(origDir);
%     ASOresult.flow.CD=0;
%     ASOresult.loop=loop;
    % ---------------------
    % Organise Outputs
    
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    obj=ASOresult.flow;
    objValue=obj.CD;
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
    
    
    % ---------------------
    % Calculate Fill Movement
    
    varExtract={'typeLoop'};
    [typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    if asoReturnFillChange
        %error('Returning Fill delta not coded yet')
        if ~isfield(ASOresult.loop,'coord')
            [ASOresult.loop.coord]=deal(ASOresult.loop.(typeLoop));
        end
        [fill,~]=LoopToFill(ASOresult.loop,baseGrid);
        additional.filldelta=fill-member.fill;
    else
        additional.filldelta=member.fill-member.fill;
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
    for ii=1:length(loop)
        
        [xMin(ii),xMax(ii),t(ii),L(ii),A(ii)]=...
            ClosedLoopProperties(loop(ii).snaxel.coord(1:end-1,:));
        
    end
    objValue=sum(A)/sum(L);
    
    additional.A=sum(A);
    additional.L=sum(L);
    additional.t=sum(t);
    additional.c=max(xMax)-min(xMin);
    additional.tc=additional.t/additional.c;
end