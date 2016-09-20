

function []=MoveToDir(nameDir,numStep)
    
     currDir=cd;
     filesepLoc=regexp(currDir,filesep);
     ii=0;
     while isempty(regexp(currDir(filesepLoc(end)+1:end),nameDir, 'once')) && ii<numStep
         
         cd ..
         currDir=cd;
         filesepLoc=regexp(currDir,filesep);
         ii=ii+1;
     end
     
    
    
    
end