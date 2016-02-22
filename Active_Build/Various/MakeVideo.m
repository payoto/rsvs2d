function []=MakeVideo(movStruct,fps,quality,datType)
    
    fileName=['Vid_',datType,'_',datestr(now,30),'.avi'];
    stampedSubFolders=['VideoArchive_',datestr(now,'yyyy_mm')];
    targetDir=[cd,'\Videos\',stampedSubFolders,'\'];
    system(['md "',targetDir,'"']);
    writerObj = VideoWriter([targetDir,fileName]);
    writerObj.FrameRate=fps;
    writerObj.Quality=quality;
    open(writerObj)
    writeVideo(writerObj,movStruct)
    close(writerObj)
end