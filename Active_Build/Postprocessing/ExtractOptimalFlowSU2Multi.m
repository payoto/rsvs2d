
function [tecPlotPre]=ExtractOptimalFlowSU2Multi(optimstruct,rootFolder,dirOptim,...
        tecPlotFile,ratio,paramoptim,allRootDir,isGradient,n)
    
    
    varExtract={'defaultVal','worker','CFDfolder'};
    [defaultVal,worker,CFDfolder]=ExtractVariables(varExtract,paramoptim);
    
    delete(tecPlotFile{1});
    delete(tecPlotFile{2});
    [iterRes,nIter,nVar]=BuildIterRes(optimstruct,defaultVal);
    
    for ii=1:numel(allRootDir)
        rootDirName{ii}=InitOptimalFlowOutput(allRootDir{ii},ratio,tecPlotFile);
    end
    rootDirName=rootDirName(~cellfun(@isempty,rootDirName));
    compType=computer;
    switch dirOptim
        case 'min'
            dir='ascend';
        case 'max'
            dir='descend';

    end
    minPos=zeros([2 0]);
    for ii=1:nIter

        itL=[optimstruct(ii).population(:).objective];
        [~,itOrd]=sort(itL,dir);
        nVarLoc=length(itL);
        iterRes(ii,1:nVarLoc)=itL;
        minPos=[minPos,[reshape(itOrd(1:min(n,nVarLoc)),[1 min(n,nVarLoc)]);...
            ones([1 min(n,nVarLoc)])]];
        %lSub1(1)=plot(ones(1,nVarLoc)*ii,iterRes(ii,1:nVarLoc),'b.','markersize',5);

    end
    
    % Prepare CFD file with newest version
    for ii=1:size(minPos,2)
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;
        %         PrepareCFDPostProcessing(minIterPos);
        
    end
    
    kk=1;
    needRerun(kk)=1;
    for ii=2:size(minPos,2)
        
        precIterPos=optimstruct(minPos(1,ii-1)).population(minPos(2,ii-1)).location;
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;
        if ~strcmp(precIterPos,minIterPos)
            kk=kk+1;
            needRerun(kk)=ii;
            
        end
        
    end
    disp([int2str(kk), ' Reruns needed, stop bitching and be patient'])
    %parfor jj=1:kk
    postList=needRerun;
    postLog=true(size(1:kk));
    for jj=1:kk
        ii=needRerun(jj);
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;

        try
            
            fileOrig=[minIterPos,filesep,'run',filesep,'flow.dat'];
            if ~exist(fileOrig,'file')
                configOrig=[minIterPos,filesep,'run',filesep,'su2.cfg'];
                configEdit=[minIterPos,filesep,'run',filesep,'su2pp.cfg'];
                fidOrig=fopen(configOrig,'r');
                fidEdit=fopen(configEdit,'w');
                while ~feof(fidOrig)
                    fprintf(fidEdit,'%s\n',regexprep(fgetl(fidOrig),...
                        '^MATH_PROBLEM=.*$','MATH_PROBLEM=DIRECT'));
                end
                fclose(fidOrig);
                fclose(fidEdit);
                system(['cd "',minIterPos,'" && SU2_SOL run',filesep,'su2pp.cfg']);
                if ~exist(fileOrig,'file')
                    error('Could not generate solution file')
                end
            end
            if isempty(FindDir([minIterPos,filesep,'run'],'flowplt_cell',false))
                
                fileNew=[minIterPos,filesep,'run',filesep,'flowplt_cell.plt'];
                fidOrig=fopen(fileOrig,'r');
                fidNew=fopen(fileNew,'w');
                str=fgetl(fidOrig);
                fprintf(fidNew,'%s \n',fgetl(fidOrig));
                fprintf(fidNew,'ZONE \n');
                fprintf(fidNew,'STRANDID = 1 \n');
                fprintf(fidNew,'SOLUTIONTIME = 1 \n');
                fprintf(fidNew,'%s \n',regexprep(fgetl(fidOrig),'ZONE ',''));
                while ~feof(fidOrig)
                    fprintf(fidNew,'%s \n',fgetl(fidOrig));
                end
                fclose(fidOrig);
                fclose(fidNew);
            end
            
        catch ME
            disp(ME.getReport);
            postLog(jj)=false;
        end
        
    end
    postList((~postLog))=[];
    minIterRootDirNum=zeros([1,size(minPos,2)]);
    
    for ii=postList
        
        minIterPos=optimstruct(minPos(1,ii)).population(minPos(2,ii)).location;
        [~,filename]=FindDir( minIterPos,'tecsubfile',false);
        jj=1;
        while isempty(regexp(minIterPos,rootDirName{jj}, 'once')) && jj<numel(rootDirName)
            jj=jj+1;
        end
        minIterRootDirNum(ii)=jj;
        jj=1;
        try
            while ~isempty(regexp(filename{jj},'copy', 'once'))
                jj=jj+1;
            end
        catch ME
            ME.getReport
            minIterPos
            optimstruct(minPos(1,ii)).population
            jj=jj
            ii=ii
            filename
            throw(ME)
        end
        
        
        copyfileRobust([minIterPos,filesep,filename{jj}],[minIterPos,filesep,...
            filename{jj},num2str(minPos(:,ii)','%i_')])
        
        copyfileRobust([[minIterPos,filesep,'run'],filesep,'flowplt_cell.plt'],...
            [[minIterPos,filesep,'run'],filesep,'flowplt_cell.plt',num2str(minPos(:,ii)','%i_')])
        
        %[snakPlt{ii}]=EditPLTTimeStrand(ii,3,2,minIterPos,[filename{jj},int2str(ii)]);
        dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'};
        expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'};
        val={(minPos(1,ii)+(minPos(2,ii)-1)/(10*(1+floor(log10(n))))),... % print decimal
            [3 4],minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5};
        nOccur=[2 2 1 1];
        [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},num2str(minPos(:,ii)','%i_')],dat,expr,val,nOccur);
        
        [flowPlt{ii}]=EditPLTTimeStrand(ii,1,2,[minIterPos,filesep,'run'],...
            ['flowplt_cell.plt',num2str(minPos(:,ii)','%i_')]);
    end
    
    % dat={'SOLUTIONTIME','STRANDID','CONNECTIVITYSHAREZONE','VARSHARELIST'}
    % expr={'SOLUTIONTIME=%f','STRANDID=%i','CONNECTIVITYSHAREZONE=%i','VARSHARELIST=([1,2]=%i)'}
    % val={ii,3,minIterRootDirNum(ii)*5,minIterRootDirNum(ii)*5}
    % nOccur=[2 2 1 1]
    % [snakPlt{ii}]=EditPLTHeader(minIterPos,[filename{jj},int2str(ii)],dat,expr,val,nOccur)
    
    for ii=postList
        
        if strcmp(compType(1:2),'PC')
            [~,~]=system(['type "',flowPlt{ii},'" >> "',tecPlotFile{1},'"']);
            [~,~]=system(['type "',snakPlt{ii},'" >> "',tecPlotFile{2},'"']);
        else
            [~,~]=system(['cat ''',flowPlt{ii},''' >> ''',tecPlotFile{1},'''']);
            [~,~]=system(['cat ''',snakPlt{ii},''' >> ''',tecPlotFile{2},'''']);
        end
    end
    tecPlotPre=regexprep(tecPlotFile,'\.plt','_pre.plt');
    if strcmp(compType(1:2),'PC')
        [d,c]=system(['preplot "',tecPlotFile{1},'" "',tecPlotPre{1},'"']);
        [d,c]=system(['preplot "',tecPlotFile{2},'" "',tecPlotPre{2},'"']);
    end
    
    
    
    
end