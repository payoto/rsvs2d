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



function []=CommitChanges(commitType,archDir2)
    
    originalDir=[cd,'\Active_Build'];
    
    timeStamp=[datestr(now,'yyyy_mm_ddHHMM')];
    stampedSubFolders=['Archive_',datestr(now,'yyyy_mm')];
    archiveDir=[cd,'\Archives\',stampedSubFolders,'\',commitType,'_',timeStamp];
    if exist('archDir2','var')
        archiveDir2=[archDir2,'\Archives\',stampedSubFolders,'\',commitType,'_',timeStamp];
        copyfile(originalDir,archiveDir2)
    end
    copyfile(originalDir,archiveDir)

end