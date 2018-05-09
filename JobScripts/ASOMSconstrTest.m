function ASOMSconstrTest(lvlSubdiv,errTreatment,startIter)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('ASOMS_subdivconstr(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter));
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/Restart_ASOMS_Short.mat'];
    ExecuteOptimisation(sprintf('ASOMS_subdivconstr(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter),{restartPath,{'none',true}});
    

end


