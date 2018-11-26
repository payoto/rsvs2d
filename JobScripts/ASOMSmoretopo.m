function ASOMSmoretopo()
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr='ASOMS_moretopo(0.2)';
    disp(callStr);
    
   
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/Restart_ASOMS_vol_0.200_moretopo.mat'];
    
    ExecuteOptimisation(callStr,{restartPath,{'none',true}});
    

end


