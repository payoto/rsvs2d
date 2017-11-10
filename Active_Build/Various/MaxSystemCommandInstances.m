function []=MaxSystemCommandInstances(launchcmd,maxExecProc,procMonit)
    % This function allows to execute system commands in DOS and
    % implementing a timeout within MATLAB after a specified amount of time
    % MATLAB continues after the call to system because a ' &' is appended
    % to the end of the command which means the command line continues the
    % execution without waiting for a result.
    %
    % INPUTS:
    %
    %   - launchcmd: The command that needs to be executed in the command
    %                line, for complex system commands use a batch file and 
    %                launch the batch file from the function
    %   - maxExecTime: The maximum execution time for the monitored process
    %                  in seconds
    %   -  procMonit: The name of process that needs to be monitored for
    %                 the timeout condition, and that needs to be stopped 
    %                 if the program doesn't complete
    
    
    launchcmd=[launchcmd,' &'];
    monitcmd=sprintf('tasklist /FI "imagename eq %s.exe" /fo table /nh',procMonit);
 
    
    [status,result] = system(launchcmd);
    %[~,~]=system(killCMDLine);
    
    % Trying to create an XFOIL timeout
    %test if the program is running
    
    [~,resultxfoil] = system(monitcmd);
    test=regexpi(resultxfoil,procMonit);
    %while the program is runnning stay in this loop
    while numel(test)>=maxExecProc
        
        [~,resultxfoil] = system(monitcmd);
        test=regexpi(resultxfoil,procMonit);

    end

    
end