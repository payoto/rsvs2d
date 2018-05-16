function ASOMSconv(lvlSubdiv,nLevel,snoptStep)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr=sprintf('ASOMS_conv(%i,%i,%i)',lvlSubdiv,nLevel,snoptStep);
    disp(callStr);
    
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/restart_ASOMS_selected.mat'];
    
    ExecuteOptimisation(callStr,{restartPath,{'none',true}});
    

end


