function [status]=CopyFileLong(pathOrig,pathNew)
    
    [origName,origPath]=ExtractName(pathOrig);
    [newName,newPath]=ExtractName(pathNew);
    if isdir(pathOrig)
        [status2,out2]=system(['robocopy "',origPath,filesep,'" "',newPath,filesep,'"']);
    else
        [status1,out1]=system(['copy "',origPath,filesep,origName,'" "',origPath,filesep,newName,'"']);
        [status2,out2]=system(['robocopy "',origPath,filesep,'" "',newPath,filesep,'" ',newName]);
        status=status2;
    end
end

function [name,path]=ExtractName(fullpath)
   
    fullpath=regexprep(fullpath,'\\\\','\\');
    posSep=regexp(fullpath,filesep);
    
    name=fullpath(posSep(end)+1:end);
    path=fullpath(1:posSep(end));
    
end