function CAFrestarts(ii)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    restartFolder = [cd,filesep,'../restarts/CAF/'];
    switch ii
        case 1
            caseStr = 'volsweeplocal(0.15,''uo'')';
            pathMat = [restartFolder,'restart_015_caf_s5.mat'];
            ExecuteOptimisation(caseStr,{pathMat, {'conjgrad',1}})
            
        case 2
            caseStr = 'volsweeplocal(0.16,''uo'')';
            pathMat = [restartFolder,'restart_016_caf_s3.mat'];
            ExecuteOptimisation(caseStr,{pathMat, {'conjgrad',1}})
        case 3
            caseStr = 'CAF_NACA0012_noref';
            ExecuteOptimisation(caseStr)
        case 4
            caseStr = 'CAF_NACA0012_ref';
            ExecuteOptimisation(caseStr)
        case 5
            caseStr = 'volsweeplocal(0.15,''uo'')';
            ExecuteOptimisation(caseStr)
            
        case 6
            caseStr = 'volsweeplocal(0.16,''uo'')';
            ExecuteOptimisation(caseStr)
    end
    
    
    
    
end