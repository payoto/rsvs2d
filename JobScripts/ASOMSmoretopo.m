function ASOMSmoretopo(vol, step)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr=sprintf('ASOMS_moretopo(%.3f,%i)',vol,step);
    disp(callStr);
    
   
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/',sprintf('Restart_ASOMS_vol_%.3f',vol),'_moretopo.mat'];
    
    ExecuteOptimisation(callStr,{restartPath,{'none',true}});
    

end


