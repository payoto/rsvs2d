
% piece of code to benchmark the flow solver options
InitialiseSnakeFlow
ExecInclude

%% Define Cases
geometriesBound=FindDir(MakePathCompliant('C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Development Notes\meshesforlaurence'),'_',1);
testCellParam(1,1:3)={'meshRefLvl',{12},1};
testCellParam(2,1:3)={'mesher',{'triangle'},1};
testCellParam(3,1:3)={'geometry',geometriesBound,0};

paramoptim=StructOptimParam('FlowSolverBench');


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
nRm=zeros(size(1,nTest));
parfor ii=1:nTest
    
      
        
         teststruct(ii).flowpath=[teststruct(ii).geometry];
%         mkdir(teststruct(ii).flowpath);
%         copyfile(teststruct(ii).geometry,teststruct(ii).flowpath)
%         fidParam=fopen([teststruct(ii).flowpath,filesep,'param.dat'],'w');
        
        boundaryLoc=FindDir(teststruct(ii).flowpath,'boundorig',0);
        if isempty(boundaryLoc)
            boundaryLoc=FindDir(teststruct(ii).flowpath,'boundary',0);
            copyfile(boundaryLoc{1},regexprep(boundaryLoc{1},'boundary','boundorig'))
        end
        loop=BoundaryInput(boundaryLoc{1});
        %copyfile(boundaryLoc{1},regexprep(boundaryLoc{1},'boundary','boundorig'))
        paramspline2=paramspline;
        nRm(ii)=0;
        for kk=1:numel(loop)
            
            if numel(loop(kk).coord)<10
                paramspline2.extrema=size(loop(kk).coord,1);
                loop(kk).coord=ResampleSpline(loop(kk).coord,paramsplinepre);
                loop(kk).coord2=loop(kk).coord;
            end
            loop(kk).coord=ResampleSpline(loop(kk).coord,paramspline2);
            [lengthParam,edgeLength]=LengthProfile(loop(kk).coord([1:end,1],:));
            loop(kk).coord(find(edgeLength(2:end)<1e-5),:)=[];
            nRm(ii)=nRm(ii)+numel(sum(edgeLength(2:end)));
        end
        
        boundaryLoc=FindDir(teststruct(ii).flowpath,'boundary',0);
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
        
        try 
            teststruct(ii).res=CutCellFlow_Handler(paramoptim1,teststruct(ii).flowpath);
%             if isempty(boundaryLoc)
%                 error('Flow Bound failed')
%             end
        catch ME
            MEflow{ii}=ME;
        end
    
    
end
nRm
