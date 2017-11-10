function BP3_RefInvDes_sweep(caseStr)
MoveToDir('source',1)
InitialiseSnakeFlow;



ExecuteOptimisation(['bp3_refsweep_',caseStr]);

end