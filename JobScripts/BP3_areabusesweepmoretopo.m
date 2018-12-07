function BP3_areabusesweepmoretopo(ii, nPop)
    if ~exist('nPop','var'); nPop = 0;end
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    callStr=sprintf('areahalfbusesweepmoretopo(%.3f,%i)', ii, nPop);
    disp(callStr);
    ExecuteOptimisation(callStr);
    



end


