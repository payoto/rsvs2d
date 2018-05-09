function [meshFile]=ASORemesh(paramoptim,dirMesh,newMesh,surfaceGeometry)
    % This function is here to support the remeshing of triangular meshes
    % using the RSVS mesh handler
    % This function takes in
    % RSVS parameter structure paramoptim
    % dirMesh : A scratch directory where data can be generated
    % newMesh : A file path where the new meshfile must be
    % surfaceGeometry: the current polygon geometry object

    % current vertices
    vertices = surfaceGeometry.vertices;

    % vertices in the order in which they need to be in the mesh
    vertConn = surfaceGeometry.getLoops();
    
    for ii=1:numel(vertConn)
        if vertConn{ii}(end)==vertConn{ii}(1)
            vertConn{ii}(end)=[];
            
        end
        loopRemesh(ii).coord=vertices(vertConn{ii},:);
        loopRemesh(ii).isccw=CCWLoop(loopRemesh.coord);
        
        loopRemesh(ii).vertorder=vertConn{ii};
    end
    % Generate triangle poly file
    
    [meshFile]=RemeshTriangle_Handler(paramoptim,dirMesh,loopRemesh);
    
    
    % Finish process
    stat=copyfile(meshFile,newMesh,'f');
    if ~stat
        error('The mesh could not be copied')
    end
    rmdir([dirMesh,filesep,'SU2CFD'],'s')
end


function []=remeshing()
    %% Set paths
    pathStr='C:\Users\ap1949\Local Documents\PhD\res\aso_1804\debug\remesh\profile_72';
    logFile=[pathStr,filesep,'log.mat'];
    workspaceFile = [pathStr, filesep, 'workspace.mat'];
    meshFile=[pathStr,filesep,'run',filesep,'mesh.su2'];
    dirMesh=[pathStr,filesep,'run'];
    newMesh=[dirMesh,filesep,'newmesh.su2'];
    %% Extract data
    
    load(logFile,'dataLog');
    load(workspaceFile,'surface');
    
    paramoptim=StructOptimParam('ASOMS_subdiv(1,''basis'',5)');
    
    vertices=dataLog.stageData(end).majorData(end-1).surface;
    [meshFile]=ASORemesh(paramoptim,dirMesh,newMesh,surface,vertices);
end
