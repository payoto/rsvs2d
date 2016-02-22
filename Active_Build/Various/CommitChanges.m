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



function []=CommitChanges(commitType1,commitType2,archDir2)
    
    originalDir=[cd,'\Active_Build'];
    
    timeStamp=[datestr(now,'yyyy_mm_ddHHMM')];
    stampedSubFolders=['Archive_',datestr(now,'yyyy_mm')];
    archiveDir=[cd,'\Archives\',stampedSubFolders,'\',commitType1,'_',timeStamp,'_',commitType2];
    if exist('archDir2','var')
        archiveDir2=[archDir2,'\Archives\',stampedSubFolders,'\',commitType1,'_',timeStamp,'_',commitType2];
        copyfile(originalDir,archiveDir2)
    end
    copyfile(originalDir,archiveDir)

end