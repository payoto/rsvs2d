function ASOMSnorestart(lvlSubdiv,errTreatment,startIter)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('ASOMS_subdiv(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter));
    restartPath='/panfs/panasas01/aero/ap1949/SnakVolParam/restarts/restart_ASOMS.mat';
    ExecuteOptimisation(sprintf('ASOMS_subdiv(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter));
    

end


