MoveToDir('source',1)
InitialiseSnakeFlow

kk=1;

try 
    ExecuteOptimisation('desk_Aero_CG_10_peak')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;
try 
    
    ExecuteOptimisation('desk_Aero_CG_10_poly')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_Aero_CG_10_polypeak')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_Aero_CG_10_none')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_SmoothCG_outmis_peak')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;
try 
    ExecuteOptimisation('desk_SmoothCG_outmis_poly')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_SmoothCG_outmis_polypeak')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;

try 
    ExecuteOptimisation('desk_SmoothCG_outmis_none')
catch ME
    errorSave{kk}=ME;
end
kk=kk+1;
