function [asoCase]=ParamConvTest(caseStr)
    
    [asoCase]=eval(['@',caseStr]);
    
end


function [options]=Default()
    
    options=ASODefaults();
    
    options.solver.np = 8;
    options.solver.mpiOpts = '';
    options.solver.mpiCmd = 'mpiexec';
    options.solver.max_iter = 1500;
    options.solver.timeout = 1800;
    options.solver.adjointMode = 'CONTINUOUS_ADJOINT';
    options.solver.mach = 2.00;
    options.solver.alpha = 0.00;
    options.solver.restartTol = 0;
    %options.solver.conf.custom. ...
    % maxFailures as parameter?
    
end

function [options]=CFLhigh()
    options=ASODefaults();
    
end

function [options]=CFLlow()
    options=ASODefaults();
    
end