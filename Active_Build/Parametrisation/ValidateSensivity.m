function [sensloops,snakloops]=ValidateSensivity(pathSens,pathSnak)
    % Takes two equivalent runs using sensitivity and snaking and returns
    % their loops
    
    
    include_Utilities;
    
    
    [sensloops]=FindProfileLoops(pathSens);
    [snakloops]=FindProfileLoops(pathSnak);
    
end

%% Loading functions
function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    subDir=dir(rootDir);
    subDir(1:2)=[];
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
end

function [profloops]=FindProfileLoops(rootDir)
    
    [profilePath,profileName]=FindDir(rootDir,'profile',true);
    profNum=regexp(profileName,'profile_','split');
    %profNum=profNum(:,2);
    
    for ii=1:length(profilePath)
        
        [restartPath,restartName]=FindDir(profilePath{ii},'restart',false);
        load(restartPath{1})
        
        profloops(str2num(profNum{ii}{2})).loop=loop;
    end
    
    
    
end

%% comparison functions

function []=CalculateDifferences(sensloops,snakloops)
    
    if numel(sensloops)~=numel(snakloops)
        error('Size mismatch beween the two sources')
    end
    
    figure,
    hold on
    
    for ii=1:numel(sensloops)
        
        errstruct(ii).coordsens=sensloops(ii).subdivspline;
        errstruct(ii).coordsnak=snakloops(ii).subdivspline;
        errstruct(ii).diff=errstruct(ii).coordsens-errstruct(ii).coordsnak;
        errstruct(ii).dist=sqrt(sum(errstruct(ii).diff.^2,2));
        errstruct(ii).sum=sum(errstruct(ii).dist);
        errstruct(ii).mean=mean(errstruct(ii).dist);
        errstruct(ii).std=std(errstruct(ii).dist);
        errstruct(ii).min=min(errstruct(ii).dist);
        errstruct(ii).max=max(errstruct(ii).dist);
        
        plot(snakloops(ii).subdivspline(:,1),errstruct(ii).dist)
    end
    
    figure,
    hold on
    plot([errstruct(:).sum])
    figure,
    hold on
    plot([errstruct(:).mean])
    figure,
    hold on
    plot([errstruct(:).std])
    figure,
    hold on
    plot([errstruct(:).min])
    figure,
    hold on
    plot([errstruct(:).max])
    
    
    
end