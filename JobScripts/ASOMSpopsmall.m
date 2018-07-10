function ASOMSpopsmall(restartPos)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr=sprintf('ASOMS_popsmall(%i)',restartPos);
    disp(callStr);
    
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/restart_ASOMS_small.mat'];
    
    ExecuteOptimisation(callStr,{restartPath,{'none',true}});
    

end


