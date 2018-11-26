function BP3_areabusesweep(vol)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    callStr=sprintf('areabusesweep(%.3f)',vol);
    disp(callStr);
    
    [~,pathHome]=system('echo -n $HOME');
    restartPath=[pathHome,'/SnakVolParam/restarts/restart_',...
        sprintf('areabuse_vol_%.3f',vol),'.mat'];
    
    
    ExecuteOptimisation(callStr,{restartPath,{'DE',true}});

end


