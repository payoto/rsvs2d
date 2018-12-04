function BP3_areabusesweepmoretopo(ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(['areahalfbusesweepmoretopo(',num2str(ii,'%e'),')'])
    ExecuteOptimisation(['areahalfbusesweepmoretopo(',num2str(ii,'%e'),')']);
    



end


