

MoveToDir('source',1)
InitialiseSnakeFlow;

try 
    ExecuteOptimisation('bp3_Aero_DE_missile_horz');
catch ME
    ME.getReport
end
try 
    ExecuteOptimisation('bp3_Aero_DE_missile');
catch ME
    ME.getReport
end
