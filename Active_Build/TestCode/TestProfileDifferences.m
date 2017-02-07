function [err1,err2]=TestProfileDifferences(optimstruct,rootDir)
    
    
    % Edit position of profiles
    for ii=1:numel(optimstruct)
        for jj=1:numel(optimstruct(ii).population)
            [optimstruct(ii).population(jj)]=TrimProfLocation(optimstruct(ii).population(jj),rootDir);
        end
    end
    
    % Identify fill errors
    for ii=1:2:(numel(optimstruct)-mod(numel(optimstruct),2));
        [jj]=1;%min([optimstruct(ii).population(:).objective])
        
        err1(ceil(ii/2),1)=1;
        err1(ceil(ii/2),2)=jj;
        err1(ceil(ii/2),3)=(optimstruct(ii).population(jj).objective-optimstruct(ii+1).population(1).objective);
        err1(ceil(ii/2),4)=sum(abs(optimstruct(ii).population(jj).fill-optimstruct(ii+1).population(1).fill));
        err1(ceil(ii/2),5:6)=CompareProfiles(optimstruct(ii).population(jj),optimstruct(ii+1).population(1));
    end
    for ii=2:2:(numel(optimstruct))-1;
        err2(ceil(ii/2),2)=1;
        err2(ceil(ii/2),1)=jj;
        [~,jj]=min([optimstruct(ii).population(:).objective]);
        err2(ceil(ii/2),3)=(optimstruct(ii).population(jj).objective-optimstruct(ii+1).population(1).objective);
        err2(ceil(ii/2),4)=sum(abs(optimstruct(ii).population(jj).fill-optimstruct(ii+1).population(1).fill));
        err2(ceil(ii/2),5:6)=CompareProfiles(optimstruct(ii).population(jj),optimstruct(ii+1).population(1));
    end
    
    
    % compare pairs of profiles
    
end


function [profile]=TrimProfLocation(profile,rootDir)
    
    
    posDir=regexp(profile.location,'iteration');
    profile.location=[rootDir,filesep,profile.location(posDir:end)];
    profile.location=MakePathCompliant(profile.location);
    
    
end

function [err]=CompareProfiles(profile1,profile2)
    
    [loop1]=LoadProfileCoord(profile1);
    [loop2]=LoadProfileCoord(profile2);
    
    if numel(loop1)~=numel(loop2)
        warning('loops have different numbers of bodies')
        err(1)=Inf;
    else
        err(1)=0;
        err(2)=0;
        totPts=0;
        for ii=1:numel(loop1)
            iMax=min(size(loop1(ii).subdivision,1),size(loop2(ii).subdivision,1));
            if size(loop1(ii).subdivision,1)~=size(loop2(ii).subdivision,1)
                disp('loops are not same sizes')
            end
            err(1)=err(1)+sum(sqrt(sum((loop1(ii).subdivision(1:iMax,:)-loop2(ii).subdivision(1:iMax,:)).^2,2)));
            totPts=totPts+iMax;
            err(2)=max(max(sqrt(sum((loop1(ii).subdivision(1:iMax,:)-loop2(ii).subdivision(1:iMax,:)).^2,2))),err(2));
        end
        err(1)=err(1)/totPts;
    end
    
end

function [loop]=LoadProfileCoord(profile)
    
    [restartFull,restartName]=FindDir(profile.location,'restart',0);
    load(restartFull{1})
    
end


function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
%     if iscell(rootDir)
%         subDir=dir(rootDir{1});
%         subDir(1:2)=[];
%         for ii=2:numel(rootDir)
%             partsubDir=dir(rootDir{ii});
%             partsubDir(1:2)=[];
%             subDir=[subDir,partsubDir];
%         end
%     else
        subDir=dir(rootDir);
        subDir(1:2)=[];
%     end
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        disp('FindDir Could not find requested item')
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end


