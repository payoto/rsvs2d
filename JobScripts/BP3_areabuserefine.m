function BP3_areabuserefine(ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    disp(['areabuserefine(',num2str(ii,'%e'),')'])
    ExecuteOptimisation(['areabuserefine(',num2str(ii,'%e'),')']);
    

end