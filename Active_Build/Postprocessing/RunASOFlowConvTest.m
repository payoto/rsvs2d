function [addstruct,population]=RunASOFlowConvTest(pathToDir,reRunDir,iter,prof,caseStr)
    % Utility to call ASOFlow from an interactive matlab session loading
    % all the necessary data
    % pathToDir <str> : Original data and run location
    % reRunDir  <str> : Location where to place the new run of the same
    %                   case
    % isRun <logical> : Provides the option to not run ASOFlow (when false) and post
    %                   treat the existing result
    
    
    if nargin<=1
        reRunDir='.';
    end
    disp('Preparing Data')
    if nargin<=4
        caseStr='Default';
    end
    [asoCase,caseStrPrint]=ParamConvTest(caseStr);
    
    [paramPath]=FindDir([pathToDir],'FinalParam',0);
    [optimPath]=FindDir([pathToDir],'OptimRes',0);
    [gridPath]=FindDir([pathToDir,filesep,...
        'iteration_0',filesep,'profile_0'],'restart',0);
    
    
    out=load(paramPath{1},'paramoptim');
    paramoptim=out.paramoptim;
    out=load(optimPath{1},'optimstruct');
    optimstruct=out.optimstruct;
    out=load(gridPath{1},'grid');
    gridBase=out.grid.base;
    
    paramoptim.obj.aso.asoCase=asoCase;
    
    % EXTRACT CORRECT PROFILES
    population=[optimstruct(iter).population];
    population=population(repmat(reshape(prof,[1 numel(prof)]),[1 numel(iter)]));
    for ii=1:numel(population)
        [loopPath,loopName]=FindDir(population(ii).location,'restart',0);
        out=load(loopPath{1},'loop');
        population(ii).loop=out.loop;
        population(ii).loopName=loopName{1};
    end
    
    % MAKE THE CORRECT RERUNDIR COPIES
    disp('Preparing directories')
    mkdir(reRunDir)
    reRunDir=[reRunDir,filesep,regexprep(pathToDir,'^.*Dir_','Dir_'),caseStrPrint];
    
    mkdir(reRunDir)
    
    for ii=1:numel(population)
        
        iterProfNum=regexprep(regexp(population(ii).loopName,'_','split'),'[^0-9]','');
        iterProfNum=cellfun(@str2double,iterProfNum(~cellfun(@isempty,iterProfNum)));
        optimstruct(iterProfNum(1)).population(iterProfNum(2));
        population(ii).oldlocation=population(ii).location;
        population(ii).location=[reRunDir,filesep,'profile_',...
            int2str(iterProfNum(1)),'_',int2str(iterProfNum(2))];
        filename=['boundary_',int2str(iterProfNum(1)),'_',int2str(iterProfNum(2)),'.dat'];
        cmd=['cp -r "',population(ii).oldlocation,filesep,filename,'" "',...
            population(ii).location,filesep,'"'];
        finDir=regexp(population(ii).location,filesep,'split');
        finDir=finDir(~cellfun(@isempty,finDir));
        population(ii).locationclean=[reRunDir,filesep,finDir{end}];
        if isdir(population(ii).location)
            system(['rm -rf ',population(ii).location])
        end
        mkdir(population(ii).location)
        system(cmd);
        disp(['Profile ',int2str(ii), ' copied'])
    end
    
    % RUN THE CONVERGENCE TESTS
    disp('Starting CFD runs for convergence testing')
    [population(:).errorMsg]=deal('');
    for ii=1:numel(population)
        origDir=cd;
        try
            [objValue,addstruct(ii)]=ASOFlowConvTest(paramoptim,population(ii),...
                population(ii).loop,gridBase);
            disp(['SUCCESS : ',population(ii).location])
        catch MEid
            disp(['FAILURE : ',population(ii).location])
            disp(MEid.getReport)
            cd(origDir)
            population(ii).constraint=false;
            population(ii).exception=[population(ii).exception,'error: ',MEid.identifier ,' - '];
            population(ii).errorMsg=[population(ii).errorMsg,MEid.getReport];
        end
    end
    save([reRunDir,filesep,'workspace_convtest.mat']);
    
    [addstruct(1:numel(population)).errorMsg]=deal(population(:).errorMsg);
    for ii=1:numel(population)
        population(ii).asocase=asoCase();
        addstruct(ii).asocase=asoCase();
    end
    save([reRunDir,filesep,'workspace_convtest.mat']);
    
    %% Additional postreatment
    
    disp('Post treatment of CFD runs')
    for jj=1:numel(population)
        try
            SU2.SU2toPLT([population(jj).locationclean,filesep,'mesh.su2'])
            SU2.SU2toPLT([population(jj).locationclean,filesep,'run',filesep,'mesh.su2'])
            
            fileToPlt={'mesh.su2.dat',['run',filesep,'mesh.su2.dat'],['run',filesep,'flow.dat']};
            
            for ii=1:numel(fileToPlt)
                cmd=['mv "',population(jj).locationclean,filesep,fileToPlt{ii},'" "',...
                    population(jj).locationclean,filesep,regexprep(fileToPlt{ii},'dat','plt'),'"'];
                system(cmd)
            end
        catch MEid2
            
        end
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
