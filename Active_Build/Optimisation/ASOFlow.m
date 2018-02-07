function [objValue,additional]=ASOFlow(paramoptim,member,loop,baseGrid)
    % Function for the interface of the RSVS shape optimisation framework
    % and the B-Spline and Smoothness constraint ASO
    %
    
    % ---------------------
    % Generate Mesh
    boundaryLoc=member.location;
    SU2Flow_Handler(paramoptim,boundaryLoc);
    copyfile([boundaryLoc,filesep,'CFDSU2',filesep,'triangularmesh.su2'],...
        [boundaryLoc,filesep,'mesh.su2'])
    rmdir([boundaryLoc,filesep,'CFDSU2'])
    
    
    % ---------------------
    % Interface parameters
    thisworker = getCurrentWorker;
    varExtract={'workerList','machineList'};
    [workerList,machineList]=ExtractVariables(varExtract,paramoptim);
    
    currentMachineFile=machineList(thisworker.ProcessId==workerList);
    optimDirectory=boundaryLoc;
    
    
    % ---------------------
    % Call ASO
    varExtract={'asoCase','asoPath','nMach'};
    [asoCase,asoPath,nMach]=ExtractVariables(varExtract,paramoptim);
    addpath(MakePathCompliant(asoPath));
    
    ASOOptions = asoCase();
    ASOOptions.solver.mach = nMach;
    ASOOptions.solver.np=currentMachineFile.slots;
    ASOOptions.solver.mpiOpts=['--hostfile ',currentMachineFile.file];
    
    ASOresult = ASO(optimDirectory,ASOOptions);
    
    
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
    if asoReturnFillChange
        %error('Returning Fill delta not coded yet')
        [fill,~]=LoopToFill(ASOresult.loop,baseGrid);
        additional.filldelta=fill-member.fill;
    else
        additional.filldelta=member.fill-member.fill;
    end
    
end