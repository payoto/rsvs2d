function BP3_TestAeroGrad(caseStr)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    steps=[1e-2,1e-3,1e-4,1e-5];
    steps=[steps;steps+1e-6;steps-1e-6];
    steps=[steps(:)]';
    for ii=steps
        ExecuteOptimisation(['TestDeriv',caseStr,'(',num2str(ii,'%.0e'),')']);
    end

end