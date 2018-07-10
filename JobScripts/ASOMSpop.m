function ASOMSpop(popNUM)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr=sprintf('ASOMS_pop(%i)',popNUM);
    disp(callStr);
    
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/restart_ASOMS.mat'];
    
    ExecuteOptimisation(callStr,{restartPath,{'none',true}});
    

end


