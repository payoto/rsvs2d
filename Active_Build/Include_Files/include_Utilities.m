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
                    (rootstruct.(fields{ii}),workingPosVec,maxPosNum);
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