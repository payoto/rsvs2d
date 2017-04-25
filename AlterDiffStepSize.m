

try
    if paramoptim.optim.CG.minDiffStep>1e-8
        paramoptim.optim.CG.diffStepSize=paramoptim.optim.CG.diffStepSize/10;
        paramoptim.optim.CG.minDiffStep=1e-9;
    end
catch MEId
    disp(MEId.getReport)
end


