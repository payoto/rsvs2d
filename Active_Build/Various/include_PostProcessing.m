function [] = include_PostProcessing()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end

%%


function []=savefig(figh,figname)

print(figh,figname,'-djpeg','-r600')

end

function [marker,t]=GenerateResultMarker(typDat)
    
    t=now;
    marker=[datestr(t,'yyyy-mm-ddTHHMMSS')...
        ,'_',typDat];
    
end

function [resultDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t)
    if ~exist('t','var'),t=now;end
    dateSubFolders=['Archive_',datestr(now,'yyyy_mm'),'\Day_',datestr(t,29)];
    resultDirectory=[resultRoot,'\',archiveName,'\',dateSubFolders,...
        '\','Dir_',marker];
    system(['md "',resultDirectory,'"']);
end