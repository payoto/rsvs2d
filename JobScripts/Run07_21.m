MoveToDir('source',1)
InitialiseSnakeFlow

kk=1;

try 
    ExecuteOptimisation('desk_Aero_CG_10_smooth_peak')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;
try 
    
    ExecuteOptimisation('desk_Aero_CG_10_smooth_poly')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_Aero_CG_10_smooth_polypeak')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_Aero_CG_10_smooth_none')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

