function BP3_AeroGrad(caseStr,ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    
    ExecuteOptimisation([caseStr,'(',num2str(ii,'%e'),')']);
    

end