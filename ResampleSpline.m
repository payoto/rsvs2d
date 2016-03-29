%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%          Spline for resampling of
%             Cloud of points
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [resampPoints]=ResampleSpline(points,paramspline)
    % Accepts a series of points and a parameter structure
    
    [parspline]=BuildSpecifiedCase(paramspline);
    
    [normPoints,parList,domSize]=GenerateParameterList(parspline,points);
    
    
    
end

%% Prepare Points

function [normPoints,parList,domSize]=GenerateParameterList(parspline,points)
    
    % Check closed curve status
    typCurve=parspline.typCurve;
    [points]=CheckClosedCurve(points,typCurve);
    
    % Normalize data
    normType=parspline.domain;
    [normPoints]=NormalizePoints(normType,points);
    
    % Extract Parameter Type
    parType=parspline.parameter;
    parList=ExtractParameterList(parType,normPoints);
    
    % Extract Domain
    [domSize]=ExtractDomain(normPoints,parList);
    
end

function [points]=CheckClosedCurve(points,typCurve)
    
   switch typCurve
       case 'closed'
           if sum(points(1,:)~=points(end,:))
               points(end+1,:)=points(1,:);
           end
       case 'open'
           
           
   end
    
    
end

function [domSize]=ExtractDomain(normPoints,parList)
    
    domSize=[min(parList),max(parList)];
    
    for ii=1:length(normPoints(1,:))
        domSize(ii+1,:)=[min(normPoints(:,ii)),max(normPoints(:,ii))];
    end
    
    
end

function [normPoints]=NormalizePoints(normType,points)
    
    switch normType
        case 'normalizeX'
            ratio=1/(max(points(:,1))-min(points(:,1)));
            
            [~,minXi]=min(points(:,1));
            vecTrans=points(minXi,:);
            
            normPoints=ratio*(points.*(ones(size(points(:,1)))*vecTrans));
            
        case 'normalizeL'
            
            [lengthParam]=LengthProfile(points);
            ratio=1/lengthParam(end);
            
            [~,minXi]=min(points(:,1));
            vecTrans=points(minXi,:);
            
            normPoints=ratio*(points.*(ones(size(points(:,1)))*vecTrans));
        otherwise
            normPoints=points;
    end
    
end

function parList=ExtractParameterList(parType,points)
    
    switch parType
        case 'x'
            parList=points(:,1);
        case 'y'
            parList=points(:,2);
        case 'l'
            parList=LengthProfile(points);
        case 'i'
            parList=(0:(length(points(:,1))-1))/(length(points(:,1))-1);
        otherwise
            warning('Unsupported Parameter for Spline Generation, using x')
            parList=points(:,1);
    end
    
    
    
end

function [lengthParam]=LengthProfile(points)
    
    points=points([1,1:end],:);
    pointsVec=points(1:end-1)-points(2:end);
    edgeLength=sqrt(sum(pointsVec.^2,2));
    lengthParam=cumsum(edgeLength);
    
    
end

%% Generate New Parameter Locations

function []=GenerateNewSampleLocations(parspline,normPoints,parList,domSize)
    
    
    newParList=Distribution(samplingDistrib,nSample,domParam);
    [extrema]=FindLocalExtremum(vec);
    
end

function [extrema]=FindLocalExtremum(vec)
    
    
    test=(~xor((vec(2:end-1)<=vec(1:end-2)),(vec(2:end-1)<=vec(3:end))))...
        *(2*((vec(2:end-1)>vec(1:end-2))-0.5));
    extrema=zeros(size(vec));
    extrema(2:end-1)=test;
    
end

function newParList=Distribution(samplingDistrib,nSample,domParam)
    
    switch samplingDistrib
        case 'even'
            newParList=linspace(domParam(1),domParam(2),nSample);
        case 'cosine'
            theta=linspace(0,pi,nSample);
            newParList=(1+cos(theta))/2*(domParam(2)-domParam(1))+domParam(1);
    end
    
end
%% Spline Parameter Handling

function [parspline]=BuildSpecifiedCase(paramspline)
    % Extracts the standard case specified and add the modifications
    % specified by paramspline
    
    [parspline]=CaseParamSpline(paramspline.splineCase);
    paramspline.structdat=GetStructureData(paramspline);
    
    fieldsInput=fieldnames(paramspline);
    
    for ii=1:length(paramspline)
        parspline.(fieldsInput{ii})=paramspline.(fieldsInput{ii});
    end
    
    
end

function [parspline]=CaseParamSpline(caseStr)
    % Main function that allows changes
    
    [parspline]=eval(['CaseSpline_',caseStr]);
    parspline.structdat=GetStructureData(parspline);
    parspline.splineCase=caseStr;
    
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

% Cases

function [parspline]=CaseSpline_default()
    
    parspline.smoothing=0;
    
    parspline.parameter='x'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='provided';
    parspline.domain='normalizeX';
    
    parspline.samplingParam='param';
    parspline.samplingN=100;
    parspline.samplingDistrib='even';
    
    parspline.LocMinDeriv='smooth';
    parspline.LocMaxDeriv='smooth';
    
    
end


