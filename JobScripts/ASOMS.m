function ASOMS(lvlSubdiv,errTreatment,startIter)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(sprintf('ASOMS_subdiv(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter));
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/restart_ASOMS.mat'];
    ExecuteOptimisation(sprintf('ASOMS_subdiv(%i,''%s'',%i)',...
        lvlSubdiv,errTreatment,startIter),{restartPath,{'none',true}});
    

end


