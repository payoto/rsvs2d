function []=TestLiftingBuse
    % piece of code to benchmark the flow solver options
    InitialiseSnakeFlow
    ExecInclude

    %% Define Cases
    % geometriesBound=FindDir(MakePathCompliant('.\Active_Build\Sample Geometries\sampleboundaries'),'boundary',0);
    testCellParam(1,1:3)={'nAlpha',{0,1, 2,3, 4,5, 6, 7, 8, 9, 10},1};
    testCellParam(2,1:3)={'area',{0.01,0.03,0.05,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.2,0.25},0};
    testCellParam(3,1:3)={'geometry',LinRes(),0};
    testCellParam(4,1:3)={'meshRefLvl',{12, 13},1};
    testCellParam(5,1:3)={'mesher',{'cutcell'},1};
    testCellParam(6,1:3)={'nMach',{2},1};

    paramoptim=StructOptimParam('FlowSolverBench');

    %% Initialise folders
    [resultRoot]=ExtractVariables({'resultRoot'},paramoptim.parametrisation);
    [marker,t]=GenerateResultMarker('FlowSolverBench');
    [writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,'ParamValidation',t);
    diaryFile=[writeDirectory,'\Latest_Diary.log'];
    diaryFile=MakePathCompliant(diaryFile);
    fidDiary=fopen(diaryFile,'w');
    fclose(fidDiary);
    diary(diaryFile);

    %% Build Test matrix
    for ii=1:size(testCellParam,1)
        nT(ii)=numel(testCellParam{ii,2});
        defStruct(1:2,ii)={testCellParam{ii,1};[]};
    end

    teststruct=repmat(struct(defStruct{:}),[1,prod(nT)]);
    nTest=prod(nT);

    for ii=1:numel(nT)
        repRate=prod([1,nT(1:ii-1)]);
        for jj=1:nTest
            teststruct(jj).(testCellParam{ii,1})=testCellParam{ii,2}...
                {mod(ceil(jj/repRate)-1,nT(ii))+1};
        end
    end
    [teststruct(:).res]=deal(struct([]));
    [teststruct(:).flowpath]=deal('');
    isParam=find([testCellParam{:,3}]);
    %% Run Test

    MEflow=cell([1,nTest]);
    MEproc=cell([1,nTest]);
    paramspline.splineCase='smoothpts';


    paramsplinepre.splineCase='presmoothpts';


    StartParallelPool(paramoptim.general.worker,4);
    
    parfor ii=1:nTest
        
         try 
            
            teststruct(ii).flowpath=[writeDirectory,filesep,'flow_',int2str(ii)];
            mkdir(teststruct(ii).flowpath);
            % copyfile(teststruct(ii).geometry,teststruct(ii).flowpath)
            teststruct(ii).geometry(teststruct(ii).area, teststruct(ii).flowpath);
            fidParam=fopen([teststruct(ii).flowpath,filesep,'param.dat'],'w');
            
            boundaryLoc=FindDir(teststruct(ii).flowpath,'boundary',0);
            loop=BoundaryInput(boundaryLoc{1});
            paramspline2=paramspline;
            for kk=1:numel(loop)
                loop(kk).coord=loop(kk).coord+(1e-16 - rand(size(loop(kk).coord))*2e-16);
                if numel(loop(kk).coord)<10
                    paramspline2.extrema=size(loop(kk).coord,1);
                    loop(kk).coord=ResampleSpline(loop(kk).coord,paramsplinepre);
                    loop(kk).coord2=loop(kk).coord;
                end
                loop(kk).coord=ResampleSpline(loop(kk).coord,paramspline2);
            end
            fidBoundary=fopen(boundaryLoc{1},'w');
            BoundaryOutput(loop,fidBoundary,'coord',0);
            fclose(fidBoundary);
            %resCell{ii}=CutCellFlow_Handler(paramoptim,boundaryLoc)
            setCell=cell(numel(isParam),1);
            kk=1;
            for jj=isParam
                setCell{kk}=teststruct(ii).(testCellParam{jj,1});
                kk=kk+1;
            end
            paramoptim1=SetVariables({testCellParam{isParam,1}},setCell,paramoptim);
            GenerateParameterFile(fidParam,paramoptim1,now,'flowtest');
            try 
                teststruct(ii).res=CutCellFlow_Handler(paramoptim1,...
                    teststruct(ii).flowpath);
    %             if isempty(boundaryLoc)
    %                 error('Flow Bound failed')
    %             end
            catch ME
                MEflow{ii}=ME;
            end
        catch ME
            
            MEproc{ii}=ME;
        end
        
    end
    save([writeDirectory,filesep,'alldata.mat']);


    %% Postprocess results
    kk=0;
    % plots with respect to refLevelCell
    listNames=cell(0);
    % for ii=1:40
    %     indStart=ii;
    %     inds=indStart:40:numel(teststruct);
    for ii=1:(numel(teststruct)/numel(testCellParam{1,2}))
        indStart=(ii-1)*numel(testCellParam{1,2})+1;
        indEnd=(ii)*numel(testCellParam{1,2});
        inds = indStart:indEnd;
        seriesName=[regexprep(regexprep(teststruct(indStart).geometry,...
            '^.*boundary_',''),'\.dat',''),...
            '_',num2str(teststruct(indStart).area,'%.2f')];
        
        jj=0;
        flagMatch=false;
        while jj<numel(listNames) && ~flagMatch
            jj=jj+1;
            flagMatch=strcmp(listNames{jj},seriesName);
        end
        
        if ~flagMatch
            listNames{end+1}=seriesName;
            h(jj+1)=figure('Name',seriesName);
        else
            figure(jj);
        end
        fieldsRes=fieldnames(teststruct(indStart).res);
        for kk=1:numel(fieldsRes)
            subplot(3,3,kk),hold on;
            title(fieldsRes{kk})
            x=[teststruct(inds).nAlpha];
            
            y=zeros(size(inds));
            for ll=1:numel(inds)
                y(ll)=[teststruct(inds(ll)).res.(fieldsRes{kk})];
            end
            plot(x,y,'-','DisplayName',[teststruct(indStart).mesher,'_',num2str(teststruct(indStart).area),'_',int2str(kk)])
        end
        if numel(fieldsRes)< 9
            subplot(3,3,numel(fieldsRes)+1),hold on;
            cl=zeros(size(inds));
            cd=zeros(size(inds));
            for ll=1:numel(inds)
                cl(ll)=[teststruct(inds(ll)).res.cl];
                cd(ll)=[teststruct(inds(ll)).res.cd];
            end
            plot(cd,cl,'-','DisplayName',[teststruct(indStart).mesher,'_',num2str(teststruct(indStart).area),'_',int2str(kk)])
            title('drag polar')
        end
    end

    for ii=1:numel(h)
        hgsave(h(ii),[writeDirectory,filesep,h(ii).Name,'.fig'])
    end

end

function [loopFunc]=LinRes()
    xMin = -0.2;
    xMax = 0.8;
    nPoints = 400;
    nMach = 2;

    loopFunc{1}=@(A, flowpath) BoundaryOutput(ConstantArea_Parabola(xMin,xMax,A,nPoints),...
        fopen([flowpath,filesep,matlab.lang.makeValidName(['boundary_Parabola',num2str(A)]),'.dat'],'w'));
%     loopFunc{2}=@(A, flowpath) BoundaryOutput(ConstantArea_Wedge(xMin,xMax,A),...
%         fopen([flowpath,filesep,matlab.lang.makeValidName(['boundary_Wedge',num2str(A)]),'.dat'],'w'));
    loopFunc{2}=@(A, flowpath) BoundaryOutput(ConstantArea_Klunker(xMin,xMax,A,nMach,nPoints),...
        fopen([flowpath,filesep,matlab.lang.makeValidName(['boundary_Klunker',num2str(A)]),'.dat'],'w'));
    loopFunc{3}=@(A, flowpath) BoundaryOutput(ConstantArea_Busemann(xMin,xMax,A,nMach),...
        fopen([flowpath,filesep,matlab.lang.makeValidName(['boundary_Busemann',num2str(A)]),'.dat'],'w'));
end