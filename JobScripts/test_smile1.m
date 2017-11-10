
MoveToDir('source',1)
InitialiseSnakeFlow;

try 
    ExecuteOptimisation('Tbp3_Aero_DE_smile_weak');
catch ME
    ME.getReport
end
try 
    ExecuteOptimisation('Tbp3_Aero_CG_smile_in');
catch ME
    ME.getReport
end
try 
    ExecuteOptimisation('Tbp3_Aero_CG_smile_out');
catch ME
    ME.getReport
end