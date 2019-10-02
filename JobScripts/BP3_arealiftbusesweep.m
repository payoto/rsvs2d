function BP3_arealiftbusesweep(vol, cl)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    if nargin>0
    	callStr=sprintf('areabusesweeplift(%.3f, %.2f)',vol,cl);
    else
    	callStr = test_areabusesweeplift;
    end
    disp(callStr);
    
    % [~,pathHome]=system('echo -n $HOME');
    % restartPath=[pathHome,'/SnakVolParam/restarts/restart_',...
    %     sprintf('areabuse_vol_%.3f',vol),'.mat'];
    
    
    ExecuteOptimisation(callStr);

end


