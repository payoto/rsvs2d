function BP3_areabuseCGRestart(ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(['areabuseCGre(',num2str(ii,'%e'),')'])
    
    
    restartDir='/panfs/panasas01/aero/ap1949/SnakVolParam/restarts/buse1705';
    restartDir='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\restarts\buse1705';
    distinct=ii;
    preStr='^OptimRes.*_areabusesweep_';
    postStr='_.mat';
    [restartPath]=IdentifyRestart(restartDir,distinct,preStr,postStr);
    
    ExecuteOptimisation(['areabuseCGre(',num2str(ii,'%e'),')'],{restartPath,{'DE',true}});
    

end