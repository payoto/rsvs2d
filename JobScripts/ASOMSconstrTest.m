function ASOMSconstrTest(lvlSubdiv,errTreatment,startIter)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('ASOMS_subdivconstr(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter));
    restartPath='/panfs/panasas01/aero/ap1949/SnakVolParam/restarts/restart_ASOMS_Short.mat';
    ExecuteOptimisation(sprintf('ASOMS_subdivconstr(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter),{restartPath,{'none',true}});
    

end


