function [] = include_Utilities()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%%

function [structdat]=ExploreStructureTree(rootstruct)
    % Recursively explores a structure to return the structures and
    % variables contained within
    
    maxPosNum=0;
    currentPosVec=[];
    [fields,vars,~]=ExploreStructureTreeRecurse(rootstruct,currentPosVec,maxPosNum);
    structdat.fields=fields;
    structdat.vars=vars;
    
    function [fields,vars,maxPosNum]=...
            ExploreStructureTreeRecurse...
            (rootstruct,currentPosVec,maxPosNum)
        % Explores a structure tree and returns
        vars=struct([]);
        fields=fieldnames(rootstruct);
        sizStruct=length(fields);
        posNumWorking=maxPosNum;
        maxPosNum=maxPosNum+sizStruct;
        
        
        for ii=1:sizStruct
            if isstruct(rootstruct.(fields{ii}))
                workingPosVec=[currentPosVec,posNumWorking+ii]; % Creates the posi
                [newfields,newvars,maxPosNum]=...
                    ExploreStructureTreeRecurse...
                    (rootstruct.(fields{ii})(1),workingPosVec,maxPosNum);
                fields=[fields;newfields];
            else
                workingPosVec=[currentPosVec,posNumWorking+ii];
                newvars.name=fields{ii};
                newvars.vec=workingPosVec;
            end
            vars=[vars,newvars];
            clear newvars
        end
        
        
    end
    
end

function [varargout]=ExtractVariables(varNames,param)
    
    varInStruct=param.structdat.vardat.names;
    varPosInStruct=param.structdat.vardat.varmatch;
    varargout{length(varNames)}=[];
    
    for ii=1:length(varNames)
        
        targLocation=regexp(varInStruct,varNames{ii}, 'once');
        if isempty(targLocation) || varPosInStruct(targLocation)==0
            error([varNames{ii},' is an invalid variable name'])
        end
        
        activeVarNum=varPosInStruct(targLocation);
        varargout{ii}=ExtractVariableValue(param,activeVarNum);
    end
    
    function [varVal]=ExtractVariableValue(param,varNum)
        
        varstruct=param.structdat.vars(varNum);
        actstruct=param;
        for jj=1:length(varstruct.vec)
            actstruct=actstruct.(param.structdat.fields{varstruct.vec(jj)});
        end
        varVal=actstruct;
    end
    
end

function [param]=SetVariables(varNames,varValues,param)
    
    varInStruct=param.structdat.vardat.names;
    varPosInStruct=param.structdat.vardat.varmatch;
    varargout{length(varNames)}=[];
    
    for ii=1:length(varNames)
        
        targLocation=regexp(varInStruct,varNames{ii}, 'once');
        if isempty(targLocation) || varPosInStruct(targLocation)==0
            error([varNames{ii},' is an invalid variable name'])
        end
        
        activeVarNum=varPosInStruct(targLocation);
        param=SetVariableValue(param,activeVarNum,varValues{ii});
    end
    
    function [param]=SetVariableValue(param,varNum,newVal)
        
        varstruct=param.structdat.vars(varNum);
        varPath='param';
        for jj=1:length(varstruct.vec)
            varPath=[varPath,'.',param.structdat.fields{varstruct.vec(jj)}];
        end
        
        eval([varPath,'=newVal;'])
    end
    
end

function [lengthParam,edgeLength]=LengthProfile(points)
    
    points=points([1,1:end],:);
    pointsVec=points(1:end-1,:)-points(2:end,:);
    edgeLength=sqrt(sum(pointsVec.^2,2));
    lengthParam=cumsum(edgeLength);
    
    
end

function [points]=RemoveIdenticalVectors(points)
    
    indOrd=1:size(points,1);
    for ii=1:size(points,2)
        [points,iRows]=SortVecColumn(points,ii);
        indOrd=indOrd(iRows);
    end
    [points,indRmv]=RemoveIdenticalConsecutivePoints(points);
    indOrd(indRmv)=[];
    [~,indOrd]=sort(indOrd);
    points=points(indOrd,:);
end

function [vec,iRows]=SortVecColumn(vec,iCol)
    % Sorts according to a columns
    [~,iRows]=sort(vec(:,iCol));
    vec=vec(iRows,:);
    
end

function [vectorEntries]=RemoveIdenticalEntries(vectorEntries)
    % Function which removes identical entries in a column vector
    vectorEntriesUnsort=vectorEntries;
    [vectorEntries,vectorIndex]=sort(vectorEntries);
    kk=1;
    rmvDI=[];
    for ii=2:length(vectorEntries)
        if vectorEntries(ii)==vectorEntries(ii-1)
            rmvDI(kk)=ii;
            kk=kk+1;
        end
    end
    %vectorEntries(rmvDI)=[];
    vectorIndex(rmvDI)=[];
    vectorEntries=vectorEntriesUnsort(vectorIndex);
end

function [trimmedPoints,indRmv]=RemoveIdenticalConsecutivePoints(points)
    
    [~,edgeLength]=LengthProfile(points);
    indRmv=find(edgeLength<1e-10);
    indRmv(1)=[]; % remove first point which has 0 distance.
    trimmedPoints=points;
    trimmedPoints(indRmv,:)=[];
    
end

function structdat=GetStructureData(paroptim)
    
    [structdat]=ExploreStructureTree(paroptim);
    structdat.vardat.names=[structdat.vars(:).name];
    structdat.vardat.varmatch=zeros(size(structdat.vardat.names));
    for ii=1:length(structdat.vars)
        jj=regexp(structdat.vardat.names,structdat.vars(ii).name);
        structdat.vardat.varmatch(jj)=ii;
    end
    
end

function CopyFileLinux(p1,p2)
    
    system(['cp -rp ''',p1,''' ''',p2,'''']);
end

function [domainBounds]=MakeCartesianGridBounds(cellLevels)
    
    %cellLevels=cellLevels+2;
    axRatio=(cellLevels+2)'/cellLevels(1);
    domainBounds=[-axRatio,axRatio];
    
end

function [domainBounds]=MakeCartesianGridBoundsInactE(cellLevels)
    
    cellNorm=cellLevels;
    cellNorm(1)=cellNorm(1)-2;
    
    cellLength=1/cellNorm(1);
    
    axRatio=(cellLevels'+2)*cellLength;
    domainBounds=[-axRatio,axRatio];
    
end

function [xMin,xMax,t,L,A]=ClosedLoopProperties(points)
    
    [A]=abs(CalculatePolyArea(points));
    vec=points([end,1:end-1],:)-points;
    L=sum(sqrt(sum(vec.^2,2)));
    t=max(points(:,2))-min(points(:,2));
    xMin=min(points(:,1));
    xMax=max(points(:,1));
    
end

function [A]=CalculatePolyArea(points)
    
    pointsVec=points';
    pointsVec=pointsVec(:);
    plot(points(:,1),points(:,2));
    n=length(points(:,1));
    centreMat=eye(2*n);
    centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;
    
    [rotDif]=[0 -1 0 1; 1 0 -1 0];
    normMat=zeros(2*n);
    for ii=1:n-1
        normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
    end
    ii=n;
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
    normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
    A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);
    
end

%{
function [A]=CalculatePolyArea(points)
    
    pointsVec=points';
    pointsVec=pointsVec(:);
    plot(points(:,1),points(:,2));
    n=length(points(:,1));
    centreMat=eye(2*n);
    centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;
    
    [rotDif]=[0 -1 0 1; 1 0 -1 0];
    normMat=zeros(2*n);
    for ii=1:n-1
        normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
    end
    ii=n;
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
    normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
    A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);
    
end
%}