function [h]=CheckOptimProfile(typeStr,varargin)
    
    
    switch typeStr
        case 'loop'
            h=figure;
            subPlotSize=[1 1];
            subPos=1;
            CheckOptimProfile_singleloopin(varargin{1},h,subPlotSize,subPos)
            
        case 'optim_all'
            rootDir=varargin{1};
            [iterPaths]=FindIter(rootDir);
            h=DisplayIterProfiles(iterPaths);
        case 'iter_all'
            iterDir=varargin{1};
            [profPaths]=FindProfile(iterDir);
            h=DisplayProfilePath(profPaths);
    end
    
    
    
    
end

function h=DisplayIterProfiles(iterPathCell)
    
    for ii=1:length(iterPathCell)
        h(ii)=DisplayProfilePath(iterPathCell{ii});
        
    end
    
end

function [h]=DisplayProfilePath(pathCell)
    
    [f1,f2]=FindCLoseFactors(length(pathCell));
    subPlotSize=[f1,f2];
    h=figure;
    for ii=1:length(pathCell)
       load(pathCell{ii});
       CheckOptimProfile_singleloopin(loop,h,subPlotSize,ii);
        
    end

end

function [iterPaths]=FindIter(rootDir)
    
    [returnPath]=FindDir(rootDir,'iteration',true);
    iterPaths={};
    for ii=1:length(returnPath)
        
        iterPaths{ii}=FindProfile(returnPath{ii});
    end

end

function [profPaths]=FindProfile(iterDir)
    
    [returnPath]=FindDir(iterDir,'profile',true);
    
    for ii=1:length(returnPath)
        
        profPaths(ii)=FindDir(returnPath{ii},'restart',false);
        
    end

end

function [f1,f2]=FindCLoseFactors(n)
    
    f=sqrt(n);
    f1=floor(f);
    f2=ceil(f);
    if f1*f2<n
        f1=f1+1;
    end
    
end

function [returnPath]=FindDir(rootDir,strDir,isTargDir)
    
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        
    end
      
    
    
end

function []=CheckOptimProfile_singleloopin(loop,h,subPlotSize,subPos)
    
    figure(h)
    subplot(subPlotSize(1),subPlotSize(2),subPos)
    hold on
    for ii=1:length(loop)
        c=loop(ii).snaxel.coord;
        plot(c(:,1),c(:,2))
        
        
        
        axis equal
        try
            c=loop(ii).subdivspline;
        catch % backwards compatibility
            c=loop(ii).subdivision;
        end
        plot(c(:,1),c(:,2))
    end
    for ii=0:2:6,
        plot([ii*1e-1 ii*1e-1],[-0.5 0.5],'k--'),
        plot(-[ii*1e-1 ii*1e-1],[-0.5 0.5],'k--'),
    end
    for ii=0:2:4,
        plot([-0.7 0.7],[ii*1e-1 ii*1e-1],'k--'),
        plot([-0.7 0.7],-[ii*1e-1 ii*1e-1],'k--'),
    end
    
end