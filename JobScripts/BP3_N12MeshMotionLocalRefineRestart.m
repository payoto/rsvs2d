function BP3_N12MeshMotionLocalRefineRestart(gridCase,optim)
    MoveToDir('source',1)
    InitialiseSnakeFlow;

    
    restartDir='/panfs/panasas01/aero/ap1949/SnakVolParam/restarts/naca1709/';
    
    
    
    distinct=[gridCase,'_',optim,'__1_2_3_4_5_6'];
    preStr='^RestartOptim.*_MMesh_';
    postStr='.mat';
    
    restartDir
    distinct
    preStr
    postStr
    [restartPath]=IdentifyRestart(restartDir,distinct,preStr,postStr);
    
    funcCall=sprintf('N12_LRef_MMesh(''%s'',''%s'')',gridCase,optim);
    disp(funcCall)
    resdat=load(restartPath);
    
    restartCell={restartPath,{optim,~mod(numel(resdat.optimstruct),2)}};
    ExecuteOptimisation(funcCall,restartCell);
    

end