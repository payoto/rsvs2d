function [objValue,additional]=ASOFlow(paramoptim,member,loop)
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
    
    currentMachineFile=machineList{thisworker.ProcessId==workerList};
    optimDirectory=boundaryLoc;
    
    
    % ---------------------
    % Call ASO
    
    
    % ---------------------
    % Organise Outputs
    
    [~,areaAdd]=LengthArea(paramoptim,member,loop);
    objValue=obj.cd;
    additional=obj;
    additional.A=areaAdd.A;
    additional.L=areaAdd.L;
    additional.t=areaAdd.t;
    additional.c=areaAdd.c;
    additional.tc=areaAdd.tc;
    
end