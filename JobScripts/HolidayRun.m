MoveToDir('source',1)
InitialiseSnakeFlow

cellError{4}=[];

try 
    
    ExecuteOptimisation('Desk_CG_missile_in');
    
catch MEexception
    
    cellError{1}=MEexception;
    
end

try 
    
    ExecuteOptimisation('Desk_CG_missile_out');
    
catch MEexception
    
    cellError{2}=MEexception;
    
end

try 
    
    ExecuteOptimisation('Desk_DE_missile_horz');
    
catch MEexception
    
    cellError{3}=MEexception;
    
end

try 
    
    ExecuteOptimisation('Desk_DE_smile_horz');
    
catch MEexception
    
    cellError{4}=MEexception;
    
end
