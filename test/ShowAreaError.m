function []=ShowAreaError(entryPoint,isfig,varargin)
    
    switch entryPoint
        
        case 'folder'
            for ii=1:size(varargin{3},1)
                [paramoptim,loop1,loop2]=FindDataFromFolder(varargin{1},...
                    varargin{2},varargin{3}(ii,:));
                [analysisCoord1,analysisCoord2{ii},h]=PlotAreaError(loop1,loop2,paramoptim);
                if ~isfig
                    
                    for jj=1:numel(h)
                        close(h(jj))
                    end
                end
            end
        case 'data'
            
            [analysisCoord1,analysisCoord2{1},h{1}]=PlotAreaError(varargin{:});
    end
    
    h=figure;
    ax(1)=subplot(3,1,1);
    hold on
    ax(2)=subplot(3,1,2);
    hold on
    ax(3)=subplot(3,1,3);
    hold on
    for ii=1:numel(analysisCoord2)
        PlotProfileDifference(analysisCoord1,analysisCoord2{ii},ax)
    end
    
end

function []=PlotProfileDifference(coord1,coord2,ax)
    
    
    if numel(coord1)==numel(coord2)
        axes(ax(1))
        deltaCoord=sqrt(sum((coord2-coord1).^2,2));
        plot(deltaCoord)
        axes(ax(2))
        deltaCoord=sqrt(((coord2-coord1).^2));
        plot(deltaCoord(:,1))
        plot(deltaCoord(:,2))
        axes(ax(3))
        plot(coord1(:,1))
        plot(coord1(:,2))
    end
    
    
end

function [analysisCoord1,analysisCoord2,h]=PlotAreaError(loop1,loop2,paramoptim)
    plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
    
    [~,h(1),~,analysisCoord1]=InverseDesign_Error(paramoptim,loop1);
    [~,h(2),targCoord,analysisCoord2]=InverseDesign_Error(paramoptim,loop2);
    
    [errorMeasure1,areaDistrib1,a1]=CompareProfilesArea(analysisCoord1,targCoord);
    [errorMeasure2,areaDistrib2,a2]=CompareProfilesArea(analysisCoord2,targCoord);
    
    h(3)=figure;
    hold on
    plotPoints(targCoord)
    plotPoints(analysisCoord1)
    plotPoints(analysisCoord2)
    
    try 
    errM=(a2(:,1)-a1(:,1))/max(abs(a2(:,1)-a1(:,1)));
    for ii=1:numel(errM);
        text(a1(ii,3),a1(ii,4),num2str(errM(ii),'%.3e'));
    end
    catch
    end
    
end

function [paramoptim,loop1,loop2]=FindDataFromFolder(pathStr,prof1,prof2)
    paramoptim=[];
    [paramPath,~]=FindDir(pathStr,'FinalParam',0);
    [p1Path,~]=FindDir(pathStr,['iteration_',int2str(prof1(1))],1);
    [p2Path,~]=FindDir(pathStr,['iteration_',int2str(prof2(1))],1);
    [p1Path,~]=FindDir(p1Path{1},['profile_',int2str(prof1(2))],1);
    [p2Path,~]=FindDir(p2Path{1},['profile_',int2str(prof2(2))],1);
    [p1Path,~]=FindDir(p1Path{1},'restart',0);
    [p2Path,~]=FindDir(p2Path{1},'restart',0);
    
    load(paramPath{1})
    load(p1Path{1},'loop')
    loop1=loop;
    load(p2Path{1},'loop')
    loop2=loop;
    
end

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    if ~isempty(subDir)
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    else
       returnSub=[]; 
    end
    
    if isempty(returnSub)
        disp('FindDir Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end
