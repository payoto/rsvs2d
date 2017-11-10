MoveToDir('source',1)
InitialiseSnakeFlow

kk=1;
try 
    
    ExecuteOptimisation('Test_polypeakCG_outmissile_Area')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('Test_peakCG_outmissile_Area')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('Test_polyCG_outmissile_Area')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_Aero_CG_10_smooth_peak')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;