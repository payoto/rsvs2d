function BP3_areabusesweep(ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(['areabusesweep(',num2str(ii,'%e'),')!tri'])
    ExecuteOptimisation(['areabusesweep(',num2str(ii,'%e'),')!tri']);
    



end


