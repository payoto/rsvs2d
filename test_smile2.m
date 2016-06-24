


InitialiseSnakeFlow;

try 
    ExecuteOptimisation('Tbp3_Aero_DE_smile_horz');
catch ME
    ME.getReport
end
try 
    ExecuteOptimisation('Tbp3_Aero_DE_smile');
catch ME
    ME.getReport
end
