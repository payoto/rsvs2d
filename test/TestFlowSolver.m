
% piece of code to benchmark the flow solver options
ExecInclude

%% Define Cases
geometriesBound=FindDir('.\Active_Build\Sample Geometries\sampleboundaries','boundary',0);
testCellParam(1,1:3)={'meshRefLvl',{10 11 12 13 },1};
testCellParam(2,1:3)={'mesher',{'cutcell','triangle'},1};
testCellParam(4,1:3)={'nMach',{0.85,2},1};
testCellParam(3,1:3)={'geometry',geometriesBound,0};

paramoptim=StructOptimParam('FlowSolverBench');

%% Initialise folders
[resultRoot]=ExtractVariables({'resultRoot'},paramoptim.parametrisation);
[marker,t]=GenerateResultMarker('FlowSolverBench');
[writeDirectory]=GenerateResultDirectoryName(marker,resultRoot,'ParamValidation',t);

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

for ii=1:nTest
    
    try 
        
        teststruct(ii).flowpath=[writeDirectory,filesep,'flow_',int2str(ii)];
        mkdir(teststruct(ii).flowpath);
        copyfile(teststruct(ii).geometry,teststruct(ii).flowpath)
        fidParam=fopen([teststruct(ii).flowpath,filesep,'param.dat'],'w');
        
        boundaryLoc=FindDir(teststruct(ii).flowpath,'boundary',0);
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
            teststruct(ii).res=CutCellFlow_Handler(paramoptim1,teststruct(ii).flowpath);
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
