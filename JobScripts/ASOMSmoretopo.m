function ASOMSmoretopo()
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr='ASOMS_moretopo';
    disp(callStr);
    
   
    
    
    ExecuteOptimisation(callStr);
    

end


