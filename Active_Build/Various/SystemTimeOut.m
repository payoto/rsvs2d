  


% execute xfoil
  procMonit='xfoil'; % Process to monitor    
  inputFile='xfoil.inp';
  outputFile='xfoil.out';
  maxExecTime=20; % in seconds
  
  wd='O:\Documents\Work\year4\research\FYP\MATLAB\Target_Functions';
  launchcmd = sprintf('cd "%s" && %s.exe < %s > %s &',wd,procMonit,inputFile,outputFile);
  monitcmd=sprintf('tasklist /FI "imagename eq %s.exe" /fo table /nh',procMonit);
  killcmd=sprintf('taskkill /im %s.exe /f',procMonit);
  
  
  killCMDLine='taskkill /im cmd.exe /f';
  
  [status,result] = system(launchcmd);
  [~,~]=system(killCMDLine);
      
  % Trying to create an XFOIL timeout
  tStart=rem(now,1)*3600*24;
  timeout=1;
  %test if the program is running
  
  [~,resultxfoil] = system(monitcmd);
  test=regexpi(resultxfoil,procMonit);
  %while the program is runnning stay in this loop
  while ~isempty(test) && timeout
        
     [~,resultxfoil] = system(monitcmd);
     test=regexpi(resultxfoil,procMonit);
      tTest=rem(now,1)*3600*24-tStart;
      % If the program takes more than  60 sec kill it 
      if tTest>maxExecTime
        [~,~]=system(killcmd);
        timeout=0;
        disp('PROCESS KILLED BY TIMEOUT')
      end
  end
  % Make sure its all closed
  [~,~]=system(killcmd);
  tTest=rem(now,1)*3600*24-tStart;
  disp(['execution completed in ',num2str(tTest,'%.4e'),' seconds'])
  %[~,~]=system(killCMDLine);