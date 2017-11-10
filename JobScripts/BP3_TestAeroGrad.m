function BP3_TestAeroGrad(caseStr)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    steps=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11];
%     steps=[steps;steps+1e-6;steps-1e-6]*10^-1;
%     steps=[steps(:)]';
    for ii=steps
        ExecuteOptimisation([caseStr,'(',num2str(ii,'%e'),')']);
    end

end