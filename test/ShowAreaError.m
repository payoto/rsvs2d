function [analysisCoord1,analysisCoord2]=ShowAreaError(entryPoint,isfig,varargin)
    modeSmoothScale='';
    switch entryPoint
        
        case 'folder'
            rootStr=varargin{1};
            iterNum=varargin{2};
            profNum=varargin{3};
            if numel(varargin)<4
                startProf=2;
            else
                startProf=varargin{4};
                
            end
            in1=[iterNum 1];
            in2=[ones([1,profNum])*iterNum;startProf:(profNum+startProf-1)]';
            for ii=1:size(in2,1)
                [paramoptim,loop1,loop2]=FindDataFromFolder(rootStr,...
                    in1,in2(ii,:));
                [analysisCoord1,analysisCoord2{ii},h]=PlotAreaError(loop1,loop2,paramoptim);
                if ~isfig
                    
                    for jj=1:numel(h)
                        close(h(jj))
                    end
                end
            end
            try
                [modeSmoothScale]=ExtractVariables({'modeSmoothScale'},paramoptim.parametrisation);
            catch
                modeSmoothScale='undef';
            end
            
            h=figure('Name',modeSmoothScale);
            pathStr=varargin{1};
            intfig=varargin{2}(1);
        case 'data'
            
            [analysisCoord1,analysisCoord2{1},h{1}]=PlotAreaError(varargin{:});
            
            h=figure;
            pathStr='./fig';
            intfig=[];
        case 'post'
            analysisCoord1=varargin{1};
            analysisCoord2=varargin{2};
            
            h=figure;
            pathStr='./fig';
            intfig=[];
    end
    
    
    ax(1)=subplot(3,1,1);
    hold on
    ylabel('Unscaled Modes')
    ax(2)=subplot(3,1,2);
    hold on
    ylabel('Normalised Modes')
    ax(3)=subplot(3,1,3);
    hold on
    ylabel('Profile')
    xlabel('Point index')
    parsplie.splineCase='inversedesign2';
    for ii=1:numel(analysisCoord2)
%         [c1,~]=ResampleSpline(analysisCoord1,parsplie);
%         [c2,~]=ResampleSpline(analysisCoord2{ii},parsplie);
        c1=analysisCoord1;
        c2=analysisCoord2{ii};
        posL=PlotProfileDifference(c1,c2,ax);
    end
    axes(ax(3))
    plot(posL,c1(:,1))
    plot(posL,c1(:,2))
    legend('x','y')
    hgsave(h,[pathStr,filesep,'iter',int2str(intfig),'_',modeSmoothScale,'.fig'])
end

function [posL]=PlotProfileDifference(coord1,coord2,ax)
    
    
    if numel(coord1)~=numel(coord2)
        
        coord1=coord1(1:min(size(coord1,1),size(coord2,1)),:);
        coord2=coord2(1:min(size(coord1,1),size(coord2,1)),:);
        
        warning('Different numbers')
    end
    axes(ax(1))
    n=size(coord1,1);
    tanVec=coord1(mod((0:n-1)+1,n)+1,:)-coord1(mod((0:n-1)-1,n)+1,:);
    posL=sqrt(sum((tanVec).^2,2));
    posL=cumsum([0;posL(1:end-1)]);
    normVec=([0 1;-1 0]*tanVec')';
    deltaCoord=sum((coord2-coord1).*normVec,2)./sqrt(sum(normVec.^2,2));
    %deltaCoord=sqrt(sum((coord2-coord1).^2,2));
    [maxVal,iMax]=max(abs(deltaCoord));
    plot(posL,deltaCoord)
    axes(ax(2))
    plot(posL,sign(iMax)*deltaCoord/maxVal)
    
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
