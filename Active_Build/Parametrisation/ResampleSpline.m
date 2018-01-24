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

function [resampPoints,splineblock]=ResampleSpline(points,paramspline)
    % Accepts a series of points and a parameter structure
    
    [parspline]=BuildSpecifiedCase(paramspline);
    points=RemoveIdenticalConsecutivePoints(points);
    points=SplitAtTrailingEdge(points,parspline.TEisLeft,parspline.samplingDistrib);
    
    [normPoints,parList,domSize]=GenerateParameterList(parspline,points);
    [nPoints,nDim]=size(normPoints);
    [splineblock]=ExtractSplineBlocks(parspline,normPoints,parList);
    [resampPoints,newPar,splineblock]=ResampleSplinePatches(splineblock,parspline,nPoints);
    [resampPoints]=ForceSpecialPoints(parspline,normPoints,parList,...
        resampPoints,newPar);
    if ExtractVariables({'splitProf'},parspline)
        resampPoints=VertSplitPointsSym(resampPoints);
    end
    
end

%% Prepare Points

function [normPoints,parList,domSize]=GenerateParameterList(parspline,points)
    
    % Check closed curve status
    typCurve=parspline.typCurve;
    [points]=CheckClosedCurve(points,typCurve);
    
    % Normalize data
    normType=parspline.domain;
    normScale=parspline.scale;
    [normPoints]=NormalizePoints(normType,points,normScale);
    
    % Extract Parameter Type
    parType=parspline.parameter;
    expectExtrema=parspline.extrema;
    parList=ExtractParameterList(parType,normPoints,expectExtrema);
    
    % Extract Domain
    [domSize]=ExtractDomain(normPoints,parList);
    
end

function [points]=CheckClosedCurve(points,typCurve)
    
    switch typCurve
        case 'closed'
            if any(points(1,:)~=points(end,:))
                points(end+1,:)=points(1,:);
            end
        case 'open'
            
        otherwise
            
            error('unknown curve type')
    end
    
    
end

function [domSize]=ExtractDomain(normPoints,parList)
    
    domSize=[min(parList),max(parList)];
    
    for ii=1:length(normPoints(1,:))
        domSize(ii+1,:)=[min(normPoints(:,ii)),max(normPoints(:,ii))];
    end
    
    
end

function [normPoints]=NormalizePoints(normType,points,normScale)
    
    switch normType
        case 'normalizeX'
            ratio=1/(max(points(:,1))-min(points(:,1)));
            
            [~,minXi]=min(points(:,1));
            vecTrans=[points(minXi,1),0];
            
            normPoints=ratio*(points-(ones(size(points(:,1)))*vecTrans));
            
        case 'normalizeL'
            
            [lengthParam]=LengthProfile(points);
            ratio=1/lengthParam(end);
            
            [~,minXi]=min(points(:,1));
            vecTrans=[points(minXi,1),0];
            
            normPoints=ratio*(points-(ones(size(points(:,1)))*vecTrans));
        case 'scaleX'
            ratio=1/normScale;
            
            [~,minXi]=min(points(:,1));
            vecTrans=[points(minXi,1),0];
            
            normPoints=ratio*(points-(ones(size(points(:,1)))*vecTrans));
            
            
        otherwise
            normPoints=points;
    end
    
end

function parList=ExtractParameterList(parType,points,expectExtrema)
    
    switch parType
        case 'x'
            parList=points(:,1);
        case 'Dx'
            parList=(abs(points(:,1)-points([1,1:end-1],1)));
            parList([false;parList(2:end)==0])=min(parList(parList~=0))/10;
            parList=cumsum(parList);
        case 'y'
            parList=points(:,2);
        case 'l'
            parList=LengthProfile(points);
        case 'clmma' % length + curvature
            if size(points,1)<6
                warning('Number of points to resample is very low and might not work with chosen distribution')
            end
            points=round(points,12);
            [parList,edgeLength]=LengthProfile(points);
            edgeLength=edgeLength/max(parList);
            parList=parList/max(parList);
            
            n=max(ceil(0.1*size(points,1)),4);
            nAverage=max(round(0.025*numel(edgeLength)),4);
            edgeCurvNorm=CurvatureParameter(points);
            
            curvParam=SmoothingParameter(parList,edgeCurvNorm,'mma',...
                nAverage,nAverage*mean(edgeLength));
            parList=parList*3/4+curvParam*1/4;
            
        case 'clint'
            if size(points,1)<6
                warning('Number of points to resample is very low and might not work with chosen distribution')
            end
            points=round(points,12);
            parList=LengthProfile(points);
            parList=parList/max(parList);
            
            edgeCurvNorm=CurvatureParameter(points);
            
            curvParam=SmoothingParameter(parList,edgeCurvNorm,'int',0,0.02);
            curvParam=(curvParam-min(curvParam))/(max(curvParam)-min(curvParam));
            parList=parList*3/4+curvParam*1/4;
        case 'clint2'
            if size(points,1)<6
                warning('Number of points to resample is very low and might not work with chosen distribution')
            end
            points=round(points,12);
            parList=LengthProfile(points);
            parList=parList/max(parList);
            
            edgeCurvNorm=CurvatureParameter(points);
            
            curvParam=SmoothingParameter(parList,edgeCurvNorm,'int2',0,0.02);
            curvParam=(curvParam-min(curvParam))/(max(curvParam)-min(curvParam));
            parList=parList*3/4+curvParam*1/4;
        case 'clintmma'
            if size(points,1)<6
                warning('Number of points to resample is very low and might not work with chosen distribution')
            end
            points=round(points,12);
            [parList,edgeLength]=LengthProfile(points);
            edgeLength=edgeLength/max(parList);
            parList=parList/max(parList);
            %             [~,edgeCurvNorm,edgeCurv,edgeCurvErr]=CurvatureProfile(points(1:end-1,:));
            %             edgeCurvNorm(end+1)=edgeCurvNorm(1);
            nAverage=max(round(0.025*numel(edgeLength)),4);
            edgeCurvNorm=CurvatureParameter(points);
            
            curvParam=SmoothingParameter(parList,edgeCurvNorm,'intmma',...
                nAverage,nAverage*mean(edgeLength));
            
            parList=parList*3/4+curvParam*1/4;
        case 'clcos'
            points=round(points,12);
            if size(points,1)<6
                warning('Number of points to resample is very low and might not work with chosen distribution')
            end
            [parList,edgeLength]=LengthProfile(points);
            edgeLength=edgeLength/max(parList);
            parList=parList/max(parList);
            
            nAverage=max(round(0.03*numel(edgeLength)),2);
            
            edgeCurvNorm=abs(LineCurvature2D(points([end-1,1:end,2],:)));
            edgeCurvNorm([1,end],:)=[];
            edgeCurvNormTest=edgeCurvNorm;
            [ptsToSave]=FindNMax(edgeCurvNormTest(1:end-1),expectExtrema);
            edgeCurvNorm(edgeCurvNorm>edgeCurvNorm(ptsToSave(end,1)))=edgeCurvNorm(ptsToSave(end,1));
            
            edgeCurvNorm=MovingAverageLoop(edgeCurvNorm',nAverage)';
            edgeCurvNorm=MovingAverageLoop(edgeCurvNorm',nAverage)';
            
            curvParam=SmoothingParameter(parList,edgeCurvNorm,'cos',expectExtrema,0.5);
%             subplot(1,2,1)
%             hold on
%             %plot(linspace(0,1,numel(curvParam)),curvParam)
%             plot(linspace(0,1,numel(curvParam)),edgeCurvNorm/max(edgeCurvNorm))
%             subplot(1,2,2)
%             hold on
%             plot(parList,curvParam)
            parList=curvParam;
        case 'i'
            parList=(0:(length(points(:,1))-1))/(length(points(:,1))-1);
        otherwise
            warning('Unsupported Parameter for Spline Generation, using x')
            parList=points(:,1);
    end
    %      assignin('base',['curvParam',int2str(kk)],edgeCurvNorm);
    %      assignin('base',['curvParamErr',int2str(kk)],edgeCurvErr);
    
end
%% Curvature parameters

function edgeCurvNorm=CurvatureParameter(points)
    edgeCurvNorm=abs(LineCurvature2D(points([end-1,1:end,2],:)));edgeCurvNorm([1,end],:)=[];
    
    edgeCurvNormTest=edgeCurvNorm;
%     edgeCurvNormTest=MovingAverageLoop(edgeCurvNormTest',5)';
%     edgeCurvNormTest=MovingAverageLoop(edgeCurvNormTest',5)';
    %figure,hold on, plot(edgeCurvNorm), plot(edgeCurvNormTest)
    
    [ptsToSave]=FindNMax(edgeCurvNormTest,4);
    edgeCurvNorm(edgeCurvNorm>edgeCurvNorm(ptsToSave(end,1)))=edgeCurvNorm(ptsToSave(end,1));
end


function [ptsToSave]=FindNMax(quantTest,nMax,tol)
    % Finds nMax local maximums in function quantTest
    % tol allows to find the remaining maximums which are at least
    % min(ptsToSave(:,2))*tol high.
    quantWeak=quantTest;
    isLocalMax=quantTest>quantTest([end,1:end-1]) &...
        quantTest>quantTest([2:end,1]);
    isLocalMaxWeak=quantTest>=quantTest([end,1:end-1]) &...
        quantTest>=quantTest([2:end,1]);

    quantTest(~isLocalMax)=-inf;
    quantWeak(~isLocalMaxWeak)=-inf;
    brkFlag=false;
    if isfinite(nMax)
        ptsToSave=zeros([nMax,2]);
        for ii=1:min(nMax,numel(quantTest))
            [ptsToSave(ii,2),ptsToSave(ii,1)]=max(quantTest);
            if ~isfinite(ptsToSave(ii,2))
                brkFlag=true;
                [ptsToSave(ii,2),ptsToSave(ii,1)]=max(quantWeak);
                if ~isfinite(ptsToSave(ii,2))
                    ptsToSave(ii:end,:)=[];
                    break
                end
            end
            quantTest(ptsToSave(ii,1))=-inf;
            quantWeak(ptsToSave(ii,1))=-inf;
        end
    else
        ptsToSave=[find(isLocalMax),quantTest(isLocalMax)];
    end
    
    if exist('tol','var')
        if brkFlag
            addPts=find(min(ptsToSave(:,2))*tol<quantWeak);
            ptsToSave=[ptsToSave;[reshape(addPts,[numel(addPts) 1]),...
                reshape(quantWeak(addPts),[numel(addPts) 1])]];
        else
            addPts=find(min(ptsToSave(:,2))*tol<quantTest);
            ptsToSave=[ptsToSave;[reshape(addPts,[numel(addPts) 1]),...
                reshape(quantTest(addPts),[numel(addPts) 1])]];
        end
    end
end

function curvParam=SmoothingParameter(parList,curvParam,type,nMMA,lInt)
    
    switch type
        case 'int2'
            
            [curvParam]=(MovingIntegralWindowLoop2(parList,(curvParam),lInt));
            [curvParam]=(MovingIntegralWindowLoop2(parList,(curvParam),lInt));
            %[curvParam]=(MovingIntegralWindowLoop2(parList,(curvParam),lInt));
        case 'int'
            
            [curvParam]=(MovingIntegralWindowLoop(parList,(curvParam),lInt));
            [curvParam]=(MovingIntegralWindowLoop(parList,(curvParam),lInt));
        case 'mma'
            curvParam=MovingAverageLoop(curvParam',nMMA)';
            curvParam=MovingAverageLoop(curvParam',nMMA)';
        case 'cos'
            [ptsToSave]=FindNMax(curvParam(1:end-1),nMMA,lInt);
            ind=sort(ptsToSave(:,1));
            nMMA=numel(ind);
            curvParam2=zeros(size(curvParam));
            
            for ii=1:numel(ind)
                indArray=(ind(ii):(ind(mod(ii,nMMA)+1)));
                
                if isempty(indArray)
                    indArray=[ind(ii):numel(parList),1:(ind(mod(ii,nMMA)+1))];
                    multiPar=(parList(end)-parList(indArray(1)))...
                        +(parList(indArray(end))-parList(1));
                    startPar=parList(indArray(1));
                    tempPar=[parList(ind(ii):numel(parList));parList(end)+parList(1:(ind(mod(ii,nMMA)+1)))];
                    tempPar=(tempPar-parList(ind(ii)))/(max(tempPar)-min(tempPar))*-2+1;
                else
                    startPar=parList(indArray(1));
                    multiPar=parList(indArray(end))-parList(indArray(1));
                    
                    tempPar=(parList(indArray)-parList(indArray(1)))/...
                        (parList(indArray(end))-parList(indArray(1)))*-2+1;
                end
%                 curvParam2(indArray)=(0.5-0.5*cos(linspace(0,pi,numel(indArray))))*...
%                     multiPar+startPar;
                %plot(linspace(0,1,numel(tempPar)),acos(tempPar)/pi)
                eps=0.98;
                if numel(tempPar)>2
                    tempPar=tempPar*eps;
                    curvParam2(indArray)=((acos(tempPar)-acos(eps))/...
                        (acos(-eps)-acos(eps)))*multiPar+startPar;
                else
                    curvParam2(indArray)=parList(indArray);
                end

            end
            curvParam=curvParam2;
            
%             curvParam(curvParam(end)<curvParam)=curvParam(curvParam(end)...
%                 <curvParam)-max(parList);%mod(curvParam2,max(parList));
             [~,b]=min(curvParam);
             if b~=1
                curvParam(1:b-1)=curvParam(1:b-1)-curvParam(b-1)+parList(b-1);
             end
            
            
        case 'intmma'
            
            [curvParam]=(MovingIntegralWindowLoop(parList,(curvParam),lInt));
            curvParam=MovingAverageLoop(curvParam',nMMA)';
            curvParam=MovingAverageLoop(curvParam',nMMA)';
            
    end
    if ~strcmp('cos',type)
        curvParam=cumsum(curvParam);
    end
    curvParam=(curvParam-min(curvParam))/(max(curvParam)-min(curvParam));
    
end

function [resampPoints]=ForceSpecialPoints(parspline,points,parList,...
        resampPoints,resampPar)
    for iii=1:size( parspline.forcePts,1)
        switch parspline.forcePts{iii,1}
            case 'maxcurv'
                edgeCurvNorm=abs(LineCurvature2D(points([end-1,1:end,2],:)));
                edgeCurvNorm([1,end],:)=[];
                edgeCurvNorm=sqrt(edgeCurvNorm)/sqrt(max(edgeCurvNorm));
                edgeCurvNormTest=edgeCurvNorm;
                [ptsToSave]=FindNMax(edgeCurvNormTest,5);
                ptsToSave(:,2)=parList(ptsToSave(:,1));
                
            case 'LETE'
                
                ptsToSave=zeros([2,2]);
                ii=1;
                [ptsToSave(ii,2),ptsToSave(ii,1)]=max(points(:,1));
                ii=2;
                [ptsToSave(ii,2),ptsToSave(ii,1)]=min(points(:,1));
                ptsToSave(:,2)=parList(ptsToSave(:,1));
            case 'all'
                ptsToSave=[(1:size(points,1))',points(:,1)];
                ptsToSave(:,2)=parList(ptsToSave(:,1));
            case 'none'
                ptsToSave=zeros([0 2]);
        end
        switch parspline.forcePts{iii,2}
            case 'add'
                for ii=1:size(ptsToSave,1)
                    ind=sum(resampPar<parList(ptsToSave(ii,1)));
                    resampPar=[resampPar(1:ind);parList(ptsToSave(ii,1));...
                        resampPar(ind+1:end)];
                    resampPoints=[resampPoints(1:ind,:);points(ptsToSave(ii,1),:);...
                        resampPoints(ind+1:end,:)];
                end
            case 'replace'
                indSave=zeros(size(ptsToSave,1),1);
                indAdd=[];
                for ii=1:size(ptsToSave,1)
                    [~,ind]=min(abs(resampPar-parList(ptsToSave(ii,1))));
                    indSave(ii)=ind;
                    if any(ind==indSave(1:ii-1))
                        indAdd=[indAdd,ii];
                    else
                        resampPar(ind)=parList(ptsToSave(ii,1));
                        resampPoints(ind,:)=points(ptsToSave(ii,1),:);
                    end
                end
                % adds points for already replaced.
                for ii=indAdd;
                    ind=sum(resampPar<parList(ptsToSave(ii,1)));
                    resampPar=[resampPar(1:ind);parList(ptsToSave(ii,1));...
                        resampPar(ind+1:end)];
                    resampPoints=[resampPoints(1:ind,:);points(ptsToSave(ii,1),:);...
                        resampPoints(ind+1:end,:)];
                end
            case 'none'
        end
    end
    [resampPoints,ptsRem]=RemoveIdenticalConsecutivePoints(resampPoints);
    resampPar(ptsRem)=[];
end

function splitPoints=SplitAtTrailingEdge(points,TEisLeft,distribution)
    splitPoints=points;
    if strcmp(distribution,'cosine') || strcmp(distribution,'2cosine') || strcmp(distribution,'split')
        TEisLeft=(TEisLeft-0.5)*2;
        [~,testGlobalMin]=min(TEisLeft*points(:,1));
        splitPoints=points([testGlobalMin:end,1:testGlobalMin-1],:);
    end
    
end

function newPoints=VertSplitPointsSym(newPoints)
    
    [xMax,ii]=max(newPoints(:,1));
    [xMin,jj]=min(newPoints(:,1));
    nP=size(newPoints,1);
    sChange=-sign(newPoints(mod(ii-1+1,nP)+1,2)-newPoints(mod(ii-1-1,nP)+1,2));
    eps=pi*1e-7/3;
    inds{1}=(min(ii,jj)+1):(max(ii,jj)-1);
    inds{2}=[1:(min(ii,jj)-1),(max(ii,jj)+1):size(newPoints,1)];
    
    indUpper=inds{1+xor(ii<jj,sChange<0)};
    indLower=inds{1+(~xor(ii<jj,sChange<0))};
    
    newPoints(indUpper,2)=newPoints(indUpper,2)+eps;
    newPoints(indLower,2)=newPoints(indLower,2)-eps;
    
    addPts=ones(2,1)*newPoints(ii,:)+[0 sChange*eps; 0 -eps*sChange];
    newPoints=[newPoints(1:ii-1,:);addPts;newPoints(ii+1:end,:)];
    
    [xMin,jj]=min(newPoints(:,1));
    
    addPts=ones(2,1)*newPoints(jj,:)+[0 -sChange*eps; 0 eps*sChange];
    newPoints=[newPoints(1:jj-1,:);addPts;newPoints(jj+1:end,:)];
    newPoints(:,1)=newPoints(:,1)-xMin+eps;
    
end

%% Separate into spline patches

function [splineblock]=ExtractSplineBlocks(parspline,normPoints,parList)
    
    nSample=parspline.samplingN;
    meanEdgeL=parspline.meanEdgeLength;
    if isempty(nSample)
        if isempty(meanEdgeL)
            nSample=size(normPoints,1);
        else % sets number of samples based on mean edge length
            [~,lProf]=LengthProfile(normPoints);
            nSample=max(ceil(sum(lProf)/meanEdgeL),21);
            
        end
    elseif nSample<0
        nSample=size(normPoints,1)*-nSample;
    end
    eps=parspline.eps;
    [interestPoints]=FindLocalExtremum(parList,eps);
    
    
    intNodes=find(sum(interestPoints,2));
    
    if ~isempty(intNodes)
        [trimBlocks,domBlocks,blockLength,totDomLength]=SplitSplineBlocks(intNodes,parList);
    else
        trimBlocks=[1,length(parList)];
        domBlocks=reshape(parList(trimBlocks),[1,2]);
        blockLength=abs(domBlocks(2)-domBlocks(1));
        totDomLength=blockLength;
    end
    
    for ii=length(trimBlocks(:,1)):-1:1
        subList=trimBlocks(ii,1):trimBlocks(ii,2);
        splineblock(ii).normPoints=normPoints(subList,:);
        splineblock(ii).parList=parList(subList);
        splineblock(ii).nNew=round(blockLength(ii)/totDomLength*nSample);
        splineblock(ii).indRange=trimBlocks(ii,:);
        splineblock(ii).domParam=domBlocks(ii,:);
    end
    
end

function [trimBlocks,domBlocks,blockLength,totDomLength]=SplitSplineBlocks(intNodes,parList)
    
    n=length(parList);
    rawBlocks=[1,intNodes(1)];
    for ii=2:length(intNodes)
        rawBlocks=[rawBlocks;[intNodes(ii-1),intNodes(ii)]];
    end
    rawBlocks=[rawBlocks;[intNodes(end),n]];
    
    nBlocks=length(rawBlocks(:,1));
    domBlocks=zeros(size(rawBlocks));
    blockLength=zeros([length(rawBlocks),1]);
    
    for ii=1:nBlocks
        domBlocks(ii,:)=parList(rawBlocks(ii,:))';
        blockLength(ii)=domBlocks(ii,2)-domBlocks(ii,1);
    end
    keepBlocks=find(blockLength~=0);
    trimBlocks=rawBlocks(keepBlocks,:);
    domBlocks=domBlocks(keepBlocks,:);
    blockLength=abs(blockLength(keepBlocks));
    totDomLength=sum(abs(blockLength));
    
    
end

function [interestPoints]=FindLocalExtremum(vec,eps)
    
    sizVec=size(vec);
    
    if sizVec(1)==1
        vec=vec';
    end
    
    
    
    testMin=(((vec(2:end-1)<=vec(1:end-2)) & (vec(2:end-1)<=vec(3:end))));
    testMax=(((vec(2:end-1)>=vec(1:end-2)) & (vec(2:end-1)>=vec(3:end))));
    testWeakSaddle=((abs(vec(2:end-1)-vec(1:end-2))<=eps) | (abs(vec(2:end-1)-vec(3:end))<=eps));
    testStrongSaddle=((abs(vec(2:end-1)-vec(1:end-2))<=eps) & (abs(vec(2:end-1)-vec(3:end))<=eps));
    
    interestPoints=zeros(length(vec),4);
    interestPoints(2:end-1,1:4)=[testMin,testMax,testWeakSaddle,testStrongSaddle];
    
end

%% Resample patches

function [newPoints,newPar,newblocks]=ResampleSplinePatches(splineblocks,parspline,nPoints)
    
    for ii=length(splineblocks):-1:1
        
        [newblocks(ii)]=ResampleBlock(splineblocks(ii),parspline,nPoints);
        
    end
    
    newPoints=vertcat(newblocks(:).newPoints);
    newPar=[newblocks(:).newParList]';
    [newPoints,indRm]=RemoveIdenticalConsecutivePoints(newPoints);
    newPar(indRm)=[];
end

function [splineblock]=ResampleBlock(splineblock,parspline,nPoints)
    %global kk
    samplingDistrib=parspline.samplingDistrib;
    nSample=splineblock.nNew;
    domParam=splineblock.domParam;
    indRange=splineblock.indRange;
    distribSource=parspline.distribution;
    provNewParList=parspline.newParList;
    typeInterp=parspline.typeInterp;
    
    splineblock.newParList=Distribution(samplingDistrib,distribSource,provNewParList,nSample,domParam);
    [cond]=PickEndcond(parspline,indRange,nPoints);
    
    [newX,splineblock.interPPX]=...
        BuildSplineInterpolant(splineblock.parList,splineblock.normPoints(:,1),...
        cond.active(1,:),cond.ValsX,splineblock.newParList,typeInterp);
    [newY,splineblock.interPPY]=...
        BuildSplineInterpolant(splineblock.parList,splineblock.normPoints(:,2),...
        cond.active(2,:),cond.ValsY,splineblock.newParList,typeInterp);
    
    splineblock.newPoints=[reshape(newX,[numel(newX) 1] ),reshape(newY,[numel(newX) 1] )];
%     assignin('base',['normPoints',int2str(kk)],splineblock.normPoints);
%     assignin('base',['newParList',int2str(kk)],splineblock.newParList);
%     assignin('base',['parList',int2str(kk)],splineblock.parList);
%     assignin('base',['newPts',int2str(kk)],splineblock.newPoints);
%     kk=kk+1;
end

function [cond]=PickEndcond(parspline,indRange,nPoints)
    % conditions are based on function csape
    % Empty condval disables the end condition and lets the spline be
    % normal
    %
    
    indepVar=parspline.parameter;
    
    switch indepVar
        case 'x'
            cond.ValsX=[0,0];
            cond.ValsY=[];
            cond.active=[1 1;0 0];
            
        case 'y'
            cond.ValsX=[];
            cond.ValsY=[0,0];
            cond.active=[0 0;1 1];
            
        case 'Dx'
            cond.ValsX=[0,0];
            cond.ValsY=[];
            cond.active=[1 1;0 0];
            
        case 'Dy'
            cond.ValsX=[];
            cond.ValsY=[0,0];
            cond.active=[0 0;1 1];
            
        otherwise
            cond.ValsX=[];
            cond.ValsY=[];
            cond.active=[1 1 ; 1 1];
            
    end
    
    if indRange(1)==1
        cond.active(:,1)=0;
    end
    if indRange(2)==nPoints
        cond.active(:,2)=0;
    end
    %cond.X2step=sum(cond.active(1,:))==1;
    %cond.Y2step=sum(cond.active(2,:))==1;
    
end

function [newPoints,interPP]=BuildSplineInterpolant(parList,data,cond,condVals,newParList,typeInterp)
    
    switch typeInterp
        case 'cubic'
            if numel(condVals)>0 && any(cond)
                if data(1)>=data(end)
                    condVals=flip(condVals);
                    cond=flip(cond);
                end
                interPP = csape(parList,[condVals(1);data;condVals(2)]',cond);
            else
                interPP = csape(parList,data');
            end
            newPoints=ppval(interPP,newParList);
        case 'linear'
            [parList,newOrd]=sort(parList);
            interPP = griddedInterpolant(parList,data(newOrd,:));
            
            
            [newPoints]=InterpolateLinear(parList,data,newParList);
        case 'shitlinear'
            
            
            [parList,newOrd]=sort(parList);
            interPP = griddedInterpolant(parList,data(newOrd,:));
            
            
            newPoints= interPP(newParList);
            
        otherwise
            error('unknown interpolation')
            
    end
    
end
function [newPoints,newPoints1]=InterpolateLinear(parList,data,newParList)
    % interpolates newPoints from newParList as defined by parList
    newPoints=zeros(size(newParList,1),size(data,2));
    %     newPoints1=zeros(size(newParList,1),size(data,2));
    nPar=numel(parList);
    ll=0;
    for ii=1:numel(newParList)
        jj=1;
        flag=false;
        while(~(flag) && jj<=nPar)
            
            flag=(parList(jj)<=newParList(ii) && parList(mod(jj,nPar)+1)>=newParList(ii));
            jjSave=jj;
            jj=jj+1;
            
        end
        if ~flag
            error('Interpolate failed as parameter is not bracketed')
        end
        
        newPoints(ii,:)=(data(mod(jjSave,nPar)+1,:)-data(jjSave,:))./...
            (parList(mod(jjSave,nPar)+1)-parList(jjSave))*(newParList(ii)-parList(jjSave))...
            +data(jjSave,:);
        %         eps=1e-5;
        %         if (abs(data(mod(jjSave,nPar)+1,:)-data(jjSave,:)) <= eps) || ...
        %              (abs((parList(mod(jjSave,nPar)+1)-parList(jjSave))) <= eps) || ...
        %              (abs((newParList(ii)-parList(jjSave))) <= eps)
        %             ll=ll+1;
        %
        %         end
        
        %         newPoints(ii,:)=(data(mod(jj,nPar)+1,:)-data(jj,:))./...
        %             (parList(mod(jj,nPar)+1)-parList(jj))*(newParList(ii)-parList(jj))...
        %             +data(jj,:);
        
    end
    
    
end

function newParList=Distribution(samplingDistrib,distribSource,provNewParList,nSample,domParam)
    
    nSample=max([nSample,2]); % so that at least the two edges are computed.
    switch distribSource
        case 'calc'
            switch samplingDistrib
                case 'even'
                    newParList=linspace(domParam(1),domParam(2),nSample);
                case 'cosine'
                    theta=linspace(0,pi,nSample);
                    newParList=(1-cos(theta))/2*(domParam(2)-domParam(1))+domParam(1);
                case '2cosine'
                    theta=linspace(0,2*pi,nSample);
                    newParList=(heaviside(theta-pi)*2-(heaviside(theta-pi)-0.5)...
                        .*2.*(1-cos(theta))/2)/2*(domParam(2)-domParam(1))+domParam(1);
                otherwise
                    error('Unknown sampling distribution')
            end
        case 'provided'
            
            newParList=provNewParList((provNewParList>=min(domParam)) & ...
                (provNewParList<=max(domParam)));
            if domParam(1)<domParam(2)
                newParList=sort(newParList,'ascend');
            else
                newParList=sort(newParList,'descend');
            end
            
    end
    
end
%% Spline Parameter Handling

function [parspline]=BuildSpecifiedCase(paramspline)
    % Extracts the standard case specified and add the modifications
    % specified by paramspline
    
    [parspline]=CaseParamSpline(paramspline.splineCase);
    parspline.structdat=GetStructureData(parspline);
    
    fieldsInput=fieldnames(paramspline);
    
    for ii=1:length(fieldsInput)
        parspline.(fieldsInput{ii})=paramspline.(fieldsInput{ii});
    end
    
    
end

function [parspline]=CaseParamSpline(caseStr)
    % Main function that allows changes
    
    [parspline]=eval(['CaseSpline_',caseStr]);
    parspline.structdat=GetStructureData(parspline);
    parspline.splineCase=caseStr;
    
end

% Cases

function [parspline]=CaseSpline_default()
    
    parspline.smoothing=0;
    parspline.TEisLeft=0;
    parspline.eps=1e-7;
    parspline.forcePts={'none','none'};
    parspline.parameter='x'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    parspline.extrema=2;
    
    parspline.distribution='provided';
    parspline.newParList=[];
    parspline.domain='normalizeX'; % 'normalizeX' 'normalizeL'
    parspline.scale=1;
    
    parspline.samplingParam='param';
    parspline.samplingN=50;
    parspline.meanEdgeLength=[];
    parspline.samplingDistrib='even';
    
    parspline.LocMinDeriv='smooth';
    parspline.LocMaxDeriv='smooth';
    parspline.typeInterp='cubic';
    parspline.splitProf=false;
    
end

function [parspline]=CaseSpline_aerofoil()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='x'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='normalizeX'; % 'normalizeX' 'normalizeL'
    
    parspline.samplingParam='param';
    parspline.samplingN=301;
    parspline.samplingDistrib='cosine';
    
end

function [parspline]=CaseSpline_snake()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='l'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='open';
    
    parspline.distribution='calc';
    parspline.domain='normalizeX'; % 'normalizeX' 'normalizeL'
    
    parspline.samplingParam='param';
    parspline.samplingN=501;
    parspline.samplingDistrib='2cosine';
    
end

function [parspline]=CaseSpline_aerosnake()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='Dx'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='normalizeX'; % 'normalizeX' 'normalizeL'
    
    parspline.samplingParam='param';
    parspline.samplingN=301;
    parspline.samplingDistrib='even';
    
end

function [parspline]=CaseSpline_inversedesign()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='x'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='normalizeX'; % 'normalizeX' 'normalizeL'
    
    parspline.samplingParam='param';
    parspline.samplingN=301;
    parspline.samplingDistrib='cosine';
    %parspline.samplingScope='local';
end

function [parspline]=CaseSpline_inversedesign2()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    parspline.eps=0;
    parspline.parameter='l'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='normalizeL'; % 'normalizeX' 'normalizeL'
    parspline.scale=1;
    
    parspline.samplingParam='param';
    parspline.samplingN=10000;
    parspline.samplingDistrib='2cosine';
    parspline.typeInterp='linear';
    %parspline.samplingScope='local';
end

function [parspline]=CaseSpline_inversedesign3()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='x'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='scaleX'; % 'normalizeX' 'normalizeL'
    parspline.scale=2.000055;
    
    parspline.samplingParam='param';
    parspline.samplingN=501;
    parspline.samplingDistrib='cosine';
    
    parspline.typeInterp='linear';
    %parspline.samplingScope='local';
end

function [parspline]=CaseSpline_naca0012()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='l'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='normalizeX'; % 'normalizeX' 'normalizeL'
    parspline.scale=1;
    
    parspline.samplingParam='param';
    parspline.samplingN=701;
    parspline.samplingDistrib='2cosine';
    parspline.splitProf=true;
    parspline.typeInterp='linear';
    %parspline.samplingScope='local';
    
end

function [parspline]=CaseSpline_convhulltri()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='l'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='none'; % 'normalizeX' 'normalizeL'
    parspline.scale=1;
    
    parspline.samplingParam='param';
    parspline.samplingN=101;
    parspline.samplingDistrib='even';
    parspline.splitProf=false;
    parspline.typeInterp='linear';
    %parspline.samplingScope='local';
    
end

function [parspline]=CaseSpline_smoothpts()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='clcos'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.forcePts={'maxcurv','replace';'LETE','add'};
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='none'; % 'normalizeX' 'normalizeL'
    parspline.scale=1;
    
    parspline.samplingParam='param';
    parspline.samplingN=[];
    parspline.meanEdgeLength=2/501;
    parspline.samplingDistrib='even';
    parspline.splitProf=false;
    parspline.typeInterp='linear';
    %parspline.samplingScope='local';
    
end
function [parspline]=CaseSpline_presmoothpts()
    
    [parspline]=CaseSpline_default();
    
    parspline.TEisLeft=0;
    
    parspline.parameter='i'; % 'y'  'l'(edge length) 'i'(index) 'Dx' (absolute change in X)
    parspline.forcePts={'all','replace';'LETE','add'};
    parspline.typCurve='closed';
    
    parspline.distribution='calc';
    parspline.domain='none'; % 'normalizeX' 'normalizeL'
    parspline.scale=1;
    
    parspline.samplingParam='param';
    parspline.samplingN=[];
    parspline.meanEdgeLength=2/51;
    parspline.samplingDistrib='even';
    parspline.splitProf=false;
    parspline.typeInterp='linear';
    %parspline.samplingScope='local';
    
end