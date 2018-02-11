function [addstruct]=RerunASOFlow(pathToDir,reRunDir,isRun)
    % Utility to call ASOFlow from an interactive matlab session loading
    % all the necessary data
    % pathToDir <str> : Original data and run location
    % reRunDir  <str> : Location where to place the new run of the same
    %                   case
    % isRun <logical> : Provides the option to not run ASOFlow (when false) and post
    %                   treat the existing result
    
    if nargin<=2
        isRun=true;
    end
    if nargin<=1
        reRunDir=pathToDir;
    end
    
    [paramPath]=FindDir([pathToDir,filesep,'..',filesep,'..'],'FinalParam',0);
    [optimPath]=FindDir([pathToDir,filesep,'..',filesep,'..'],'OptimRes',0);
    [gridPath]=FindDir([pathToDir,filesep,'..',filesep,'..',filesep,...
        'iteration_0',filesep,'profile_0'],'restart',0);
    [loopPath,loopName]=FindDir(pathToDir,'restart',0);
    
    out=load(paramPath{1},'paramoptim');
    paramoptim=out.paramoptim;
    out=load(optimPath{1},'optimstruct');
    optimstruct=out.optimstruct;
    out=load(gridPath{1},'grid');
    gridBase=out.grid.base;
    out=load(loopPath{1},'loop');
    loop=out.loop;
    
    iterProfNum=regexprep(regexp(loopName{1},'_','split'),'[^0-9]','');
    iterProfNum=cellfun(@str2double,iterProfNum(~cellfun(@isempty,iterProfNum)));
    member=optimstruct(iterProfNum(1)).population(iterProfNum(2));
    if ~strcmp(reRunDir,pathToDir)
        mkdir(reRunDir)
        cmd=['cp -r "',pathToDir,filesep,'" "',reRunDir,filesep,'"'];
        finDir=regexp(pathToDir,filesep,'split');
        finDir=finDir(~cellfun(@isempty,finDir));
        reRunDir=[reRunDir,filesep,finDir{end}];
        member.location=reRunDir;
        system(cmd);
    end
    try
        if isRun
            [objValue,addstruct]=ASOFlow(paramoptim,member,...
                loop,gridBase);
        else
            objValue=[]
            addstruct=[];
        end
    catch MEid
        
    end
    %% Additional postreatment
    try
        SU2.SU2toPLT([reRunDir,filesep,'mesh.su2'])
        SU2.SU2toPLT([reRunDir,filesep,'run',filesep,'mesh.su2'])
        
        fileToPlt={'mesh.su2.dat',['run',filesep,'mesh.su2.dat'],['run',filesep,'flow.dat']};
        
        for ii=1:numel(fileToPlt)
            cmd=['mv "',reRunDir,filesep,fileToPlt{ii},'" "',...
                reRunDir,filesep,regexprep(fileToPlt{ii},'dat','plt'),'"'];
            system(cmd)
        end
    catch MEid2
        
    end
    
    try
        load([reRunDir,filesep,'log.mat'])
        kk=0;
        for ii=2:size(dataLog,1)
            if dataLog.majorIt(ii)~=dataLog.majorIt(ii-1)
                kk=kk+1;
                h(kk)=figure;
                hold on;
                plot(dataLog.surface{ii-1}(:,1),dataLog.surface{ii-1}(:,2),'k*')
            end
            plot(dataLog.surface{ii}(:,1),dataLog.surface{ii}(:,2),'+')
        end
        for ii=1:numel(h)
            hgsave(h(ii),[reRunDir,filesep,'majoritmvmt_',int2str(ii),'.fig'])
        end
    catch
        
    end
    
    if exist('MEid','var')
        try
            OutputErrorReport(MEid,reRunDir,{pathToDir,reRunDir});
        catch
        end
        rethrow(MEid)
        
    elseif exist('MEid2','var')
        rethrow(MEid2)
    else
        disp('exited without errors')
    end
end

function []=OutputErrorReport(ME,pathOut,errPaths)
    t=now;
    filename=['bughpc_',datestr(t,'yymmdd'),'.txt'];
    
    fid=fopen([pathOut,filesep,filename],'w');
    fprintf(fid,'date : %s \n',datestr(t));
    for ii=1:numel(errPaths)
        fprintf(fid,'Path: %s \n',errPaths{ii});
    end
    fprintf(fid,'Identifier: %s \n\n',ME.identifier);
    fprintf(fid,'Frequency: \n\n Observations: \n\n Solution: \n\n');
    
    
    fprintf(fid,'Full Error Message:\n\n %s ',ME.getReport);
    fclose(fid);
    
end
