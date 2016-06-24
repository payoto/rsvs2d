


InitialiseSnakeFlow;

try 
    ExecuteOptimisation('Tbp3_Aero_DE_missile');
catch ME
    ME.getReport
end
try 
    ExecuteOptimisation('Tbp3_Aero_DE_missile_horz');
catch ME
    ME.getReport
end
