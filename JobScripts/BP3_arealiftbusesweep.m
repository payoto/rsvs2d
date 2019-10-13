function BP3_arealiftbusesweep(vol, cl, isRestart)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    if nargin < 3
        isRestart = false;
    end
        
    if nargin>0
    	callStr=sprintf('areabusesweeplift(%.3f,%.2f)',vol,cl);
    else
    	callStr = 'test_areabusesweeplift';
    end
    disp(callStr);
    
    if ~isRestart

        ExecuteOptimisation(callStr);

    else

        [~,pathHome]=system('echo -n $HOME');
        restartPath=[pathHome,'/SnakVolParam/restarts/buselift/restart_',...
            CallString2Name(callStr, 'legacy'),'.mat'];
        disp(restartPath);

        
        ExecuteOptimisation(callStr, {restartPath, {'DE', false}});
    end
end


