%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Commit Active Build to 
%                 Archive
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function []=CommitChanges(commitType)
    
    originalDir=[cd,'\Active_Build'];
    
    timeStamp=[datestr(now,'yyyy_mm_ddHHMM')];
    stampedSubFolders=['Archive_',datestr(now,'yyyy_mm')];
    archiveDir=[cd,'\Archives\',stampedSubFolders,'\',commitType,'_',timeStamp];
    
    copyfile(originalDir,archiveDir)

end