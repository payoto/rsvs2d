%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          snake parametrisation
%      for Aerodynamic shape parametrisation
%           - Outputs Management Function - 
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function []=ManageOutputResults(param,loop,tecoutstruct)
    
    % Unpack necessary variables
    varExtract={'makeMov','typDat','resultRoot','archiveName'};
    [makeMov,typDat,resultRoot,archiveName]=ExtractVariables(varExtract,param);
    
    
    % Create Marker
    [marker]=GenerateResultMarker(typDat);
    % Create Directory
    [writeDirectory]=GenerateDirectoryName(marker,resultRoot,archiveName);
    
    % Output boundary data file
    [fidBoundary]=OpenBoundaryFile(writeDirectory,marker);
    BoundaryOutput(loop,fidBoundary);
    fclose(fidBoundary);
    
    % TecPlot Data
    [fidTecPLT]=OpenTecPLTFile(writeDirectory,marker);
    fidTecLAY=OpenTecLayFile(writeDirectory,marker);
    
    baseGrid=tecoutstruct.baseGrid;
    fineGrid=tecoutstruct.fineGrid;
    snakSave=tecoutstruct.snakSave;
    connectstructinfo=tecoutstruct.connectstructinfo;
    
    TecplotOutput(fidTecPLT,baseGrid,fineGrid,snakSave,connectstructinfo)
    
    % Parameter Data
    [fidParam]=OpenParamFile(writeDirectory,marker);
    
    
    
    % Comments file
    [fidComments]=OpenCommentsFile(writeDirectory,marker);
    
    
    % Video File
    
    if makeMov
        fps=5;
        quality=100;
        movStruct=[snakSave(:).movFrame];
        [fileName]=NameVideoFile(writeDirectory,marker);
        MakeVideo(movStruct,fps,quality,fileName);
    end
    
    
    
    
    
    
    
    
    fclose('all');
    
end


%% Boundary Output to .dat file

function []=BoundaryOutput(loop,FID)
    
    % trim loops and extract data
    loopout=TrimLoops(loop);
    % format numeric data to printable strings
    cellLoops=DataToString(loopout);
    
    % print string data to file
    WriteToFile(cellLoops,FID)
    
end

function [loopout]=TrimLoops(loop)
    % function extracting the data that must be written to the boundary.dat
    % file
    
    nLoop=length(loop);
    
    for ii=1:nLoop
        if loop(ii).isccw
            loopout.surf(ii).coord=loop(ii).subdivision(1:end-2,:);
        else
            loopout.surf(ii).coord=loop(ii).subdivision(end-2:-1:1,:);
        end
        loopout.surf(ii).nvertex=length(loopout.surf(ii).coord(:,1));
        loopout.surf(ii).nfaces=loopout.surf(ii).nvertex;
    end
    loopout.total.nvertex=sum([loopout.surf(:).nvertex]);
    loopout.total.nfaces=sum([loopout.surf(:).nfaces]);
    loopout.total.nloops=nLoop;
    
end

function cellLoops=DataToString(loopout)
    % Transforms the data organised in loopout into a cell array of strings
    % ready to be written to a file
    
    cellLoops{1}=int2str(loopout.total.nloops);
    cellLoops{2}=int2str([loopout.total.nvertex, loopout.total.nfaces]);
    kk=3;
    
    % run through the different surfaces
    for ii=1:loopout.total.nloops
        cellLoops{kk}=int2str(loopout.surf(ii).nvertex);
        kk=kk+1;
        % Write an array of numbers in the correct format
        for jj=1:length(loopout.surf(ii).coord(:,1))
            cellLoops{kk}='';
            for ll=1:length(loopout.surf(ii).coord(jj,:))
                cellLoops{kk}=[cellLoops{kk},...
                    num2str(loopout.surf(ii).coord(jj,ll),8),'  '];
            end
            kk=kk+1;
        end
        
    end
    
end

function []=WriteToFile(cellLoops,FID)
    % writes the data to the file
    
    
    for ii=1:length(cellLoops)
        for jj=1:length(cellLoops{ii}(:,1))
            fprintf(FID,[cellLoops{ii}(jj,:),'\n']);
        end
    end
end

%% Video Output functions

function []=MakeVideo(movStruct,fps,quality,fileName)
    
    writerObj = VideoWriter(fileName);
    writerObj.FrameRate=fps;
    writerObj.Quality=quality;
    open(writerObj)
    writeVideo(writerObj,movStruct)
    close(writerObj)
end

%% File Opening Functions

function [FID]=OpenBoundaryFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['boundary_',marker,'.dat'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=OpenTecPLTFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['tec360dat_',marker,'.plt'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=OpenTecLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['tec360lay_',marker,'.lay'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=OpenParamFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['param_',marker,'.dat'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=OpenCommentsFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['Comments_',marker,'.txt'];
    FID=fopen([writeDirectory,'\',fileName],'w+');
    
end

function [FID]=NameVideoFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['Video_',marker,'.avi'];
    FID=[writeDirectory,'\',fileName];
    
end

%% Open Result Directory


function [marker]=GenerateResultMarker(typDat)
    
    
    marker=[datestr(now,'yyyy-mm-ddTHHMMSS')...
        ,'_',typDat];
     
end


function [resultDirectory]=GenerateDirectoryName(marker,resultRoot,archiveName)
    t=now;
    dateSubFolders=['Archive_',datestr(now,'yyyy_mm'),'\Day_',datestr(t,29)];
    resultDirectory=[resultRoot,'\',archiveName,'\',dateSubFolders,...
        '\','Dir_',marker];
    system(['md "',resultDirectory,'"']);
end





