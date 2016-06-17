function taskStartup(task)
    % TASKSTARTUP Perform user-specific task startup actions on a worker
    %
    %   taskStartup(task)
    %
    %   To define specific task initialization actions on each worker for each
    %   task, provide a taskStartup function that is accessible to the worker.
    %   This function can be in any of the following locations, searched for in
    %   this order:
    %
    %   1. In a file that is specified in the job object's FileDependencies
    %      property.
    %   2. In a folder that is specified in the job object's PathDependencies
    %      property.
    %   3. In this file.
    %
    %   The task parameter that is passed to this function is the task object
    %   that the worker is executing.
    %
    %   If this function throws an error, the error information appears in the
    %   task object's Error, ErrorMessage, and ErrorIdentifier properties, and
    %   the task function and taskFinish will not be evaluated. Any callback
    %   functions for the task will still be evaluated.
    %
    %   If the same MATLAB worker session is used for running multiple tasks or
    %   jobs, any path changes made here or during the execution of tasks will
    %   be reverted by the MATLAB Distributed Computing Server to their
    %   original values before the next job runs, but preserved for subsequent
    %   tasks in the same job. Any data stored by this function or by the
    %   execution of the job's tasks (for example, in the base workspace or in
    %   global or persistent variables) will not be cleared by the MATLAB
    %   Distributed Computing Server before the next job runs, unless the
    %   RestartWorker property of the next job is set to true.
    %
    %   See also jobStartup, taskFinish.
    
    % Copyright 2004-2012 The MathWorks, Inc.
    
    %#ok<*INUSD> Suppression because this example does not use the input argument
    
     InitialiseSnakeFlow
     include_SnakeParam
     include_EdgeInformation
     include_Utilities
     include_PostProcessing
     include_Optimisation
%     include_Mex_Wrapper