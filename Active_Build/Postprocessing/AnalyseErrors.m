function [errstruct]=AnalyseErrors(rootDir)
    %Extracts information about all the Error Report files in a directory
    
    [childFolder]=ExploreFolderTree(rootDir);
    
%    FindDir([minIterPos,filesep,'CFD'],'Tec360plt_',false)
    [pltPaths]=FindErrorFiles(childFolder);
    errstruct=repmat(struct('id','','iter',[],'prof',[],'path','','indstr','','message','','pattern',''),[1 0]);
    disp('starting file operations')
    for ii=1:numel(pltPaths)
        [errstruct]=ReadInErrors(pltPaths{ii},errstruct);
    end
    save([rootDir{1},filesep,'ErrorFilesSummary.mat'],'errstruct');
end


function [errstruct]=ReadInErrors(errPath,errstruct)
    
    fid=fopen(errPath,'r');
    currMode=1;
    while ~feof(fid)
        
        str=fgets(fid);
        startRead=~isempty(regexp(str(1:min(numel(str),1)),'[0-9]', 'once'));
        falseRead=~isempty(regexp(str(1:min(numel(str),1)),'[#,-]', 'once'));
        keeprun=true;
        while startRead || keeprun
            switch currMode
                case 1 % looks for start of error message
                    
                    if startRead
                        cellStr=regexp(str,',','split');
                        errorLoc=cellfun(@isempty,regexp(cellStr,'error'));
                        if any(errorLoc)
                            currerr.id=cellStr{find(~errorLoc,1,'first')};
                            currerr.iter=str2double(cellStr{1});
                            currerr.prof=str2double(cellStr{2});
                            currerr.path=errPath;
                            currerr.indstr=str;
                            currerr.message='';
                            currerr.pattern='';
                            currMode=2;
                            kk=0;
                        end
                    end
                    keeprun=false;
                    startRead=false;
                case 2 % Case to read the text;
                    if ~startRead && ~falseRead
                        kk=kk+1;
                        currerr.message=char(currerr.message,str);
                        if kk<=4
                            currerr.pattern=char(currerr.pattern,str);
                        end
                        keeprun=false;
                    else
                        keeprun=true;
                        currMode=3;
                    end
                case 3 % Finish the treatment of the error
                    
                    isSame=false(size(errstruct));
                    for ii=1:numel(errstruct)
                        if strcmp(errstruct(ii).id,currerr.id)
                            if strcmp(errstruct(ii).pattern,currerr.pattern)
                                isSame(ii)=true;
                                break
                            end
                        end
                    end
                    if any(isSame)
                        errstruct(isSame).iter=[errstruct(isSame).iter,currerr.iter];
                        errstruct(isSame).prof=[errstruct(isSame).prof,currerr.prof];
                        errstruct(isSame).path=char(errstruct(isSame).path,currerr.path);
                        
                        
                    else
                        
                        errstruct=[errstruct,currerr]; %#ok<AGROW>
                        
                    end
                    keeprun=false;
                    currMode=1;
            end
            
        end
                
           
        
        
    end
    
    
    
    fclose(fid);
    
end

function [pltPaths]=FindErrorFiles(childFolder)
    kk=1;
    for ii=1:length(childFolder)
        intermPath=FindDir(childFolder{ii},'ErrorReport',false);
        if ~isempty(intermPath)
            pltPaths(kk:kk+numel(intermPath)-1)=intermPath;
            kk=kk+numel(intermPath);
        end
    end
    
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
        %fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end

function [addFolders]=ExploreFolderTree(rootDir)
    % adds a set of paths to the active path
    addFolders=rootDir;
    for ii=1:length(rootDir)
        dirinfo=dir(rootDir{ii});
        dirNames={dirinfo([dirinfo(:).isdir]).name};
        dirNames(1:2)=[];
        if numel(dirNames)>0
            branchDir={''};
            for jj=1:length(dirNames)
                branchDir{jj}=[rootDir{ii},filesep,dirNames{jj}];
            end
            
            [addSubFolders]=ExploreFolderTree(branchDir);
            addFolders=[addFolders,addSubFolders];
        end
    end
    
end
