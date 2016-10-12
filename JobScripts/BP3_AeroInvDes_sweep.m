function BP3_AeroInvDes_sweep(caseStr)
MoveToDir('source',1)
InitialiseSnakeFlow;



ExecuteOptimisation(caseStr);

end