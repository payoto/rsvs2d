function ASOMSvol(vol)
    
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    callStr=sprintf('ASOMS_vol(%.3f)',vol);
    disp(callStr);
    
    
    
    
    ExecuteOptimisation(callStr);
    

end


