function ASOMSvol(vol)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr=sprintf('ASOMS_vol(%.3f)',vol);
    disp(callStr);
    
    
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/restart_',...
        sprintf('ASOMS_vol_%.3f',vol),'.mat'];
    
    
    ExecuteOptimisation(callStr,{restartPath,{'none',true}});
    

end


