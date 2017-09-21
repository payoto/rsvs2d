
MoveToDir('source',1)
InitialiseSnakeFlow;

try 
    ExecuteOptimisation('bp3_Aero_CG_smile_in');
catch ME
    ME.getReport
end
try 
    ExecuteOptimisation('bp3_Aero_CG_smile_out');
catch ME
    ME.getReport
end