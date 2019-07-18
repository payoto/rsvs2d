%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%            - Grid Initialisation -
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [unstructured,loop,unstructReshape]=...
        GridInitialisationV2(param)
    % Main function for the execution of the Subdivision process
    
    %unpacking input parameters
    varExtract={'passDomBounds','passGridSteps','passPadding',...
        'typDat','loadLogical','isCheckRes','boundstr','gridDistrib'};
    [passDomBounds,passGridSteps,passPadding,...
        typDat,loadLogical,isCheckRes,boundstr,gridDistrib]=ExtractVariables(varExtract,param);
    
    % Defining global variables
    
    global domainBounds % domainSize is the Bounds of the design space Form [Xmin,Xmax;Ymin,Ymax;...]
    domainBounds=passDomBounds;
    
    global nGridSteps % number of steps in design domain
    nGridSteps=passGridSteps;
    
    global nDim % number of dimensions
    nDim=length(domainBounds(:,1));
    
    
    global nPadding
    nPadding=passPadding;
    
    %     typDat='foillo';
    %     loadLogical=false;
    %     isCheckRes=true;
    
    % Execution of subroutines
    if ~loadLogical
        
        [unstructured]=InitialisEdgeGrid(param);
        [unstructured]=GridRedistrib(unstructured,gridDistrib);
        [unstructured]=EdgeProperties(unstructured);
        isEdge=unstructured.edge.(boundstr{1});
        cond=boundstr{3};
        [loop]=OrderSurfaceVertex(unstructured,isEdge,cond);
    else
        datPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.mat'];
        load(datPath);
        if ~exist('loop','var')
            isEdge=unstructured.edge.(boundstr{1});
            cond=boundstr{3};
            [loop]=OrderSurfaceVertex(unstructured,isEdge,cond);
        end
    end
    
    if isCheckRes
        %CheckResults(unstructured,loop)
    end
    disp('Reshaping')
    [unstructReshape]=ModifUnstructured(unstructured);
    [loop]=EdgeInCondForVertex(loop,unstructReshape,cond);
    [unstructReshape] = FillNewGrid(unstructReshape, param);
    sub = FindObjNum([],unstructured.cell.index,[unstructReshape.cell.index]);
    unstructured.cell.fill = [unstructReshape.cell(sub).fill];
    
end

function [unstructured]=GridRedistrib(unstructured,gridDistrib)
    
    switch gridDistrib
        
        case 'none'
            
        case 'cosX1'
            xMax=1; %max(coord(:,1));
            xMin=-1; %min(coord(:,1));
            [unstructured]=CosGridDistribX(unstructured,xMax,xMin);
        case 'cosX01'
            xMax=1; %max(coord(:,1));
            xMin=0; %min(coord(:,1));
            [unstructured]=CosGridDistribX(unstructured,xMax,xMin);
        case 'cosXsquared01'
            xMax=1; %max(coord(:,1));
            xMin=0; %min(coord(:,1));
            [unstructured]=Cos2GridDistribX(unstructured,xMax,xMin);
        case 'cosXYsquared01'
            xMax=1; %max(coord(:,1));
            xMin=0; %min(coord(:,1));
            [unstructured]=Cos2GridDistribX(unstructured,xMax,xMin);
            [unstructured]=CosGridDistribY(unstructured);
        case 'cosXsquared'
            xMax=max(unstructured.vertex.coord(:,1));
            xMin=min(unstructured.vertex.coord(:,1));
            [unstructured]=Cos2GridDistribX(unstructured,xMax,xMin);
        case 'cosX'
            xMax=max(unstructured.vertex.coord(:,1));
            xMin=min(unstructured.vertex.coord(:,1));
            [unstructured]=CosGridDistribX(unstructured,xMax,xMin);
        case 'bilinearcos'
            tol = 0.0002; % tolerance 
            centre=0.35;
            xMax=1-tol/2;
            xMin=0+tol/2;
            [unstructured]=HalfCosGridDistribX(unstructured,xMax,xMin);
            [unstructured]=BiLinearGridDistribX(unstructured,xMax,xMin,centre);
        case 'bilinear'
            centre=0.3;
            minX=0;
            maxX=1;
            [unstructured]=BiLinearGridDistribX(unstructured,maxX,minX,centre);
        case 'thinLETE'
            maxExcess=0.1;
            xMax=1;
            xMin=0;
            [unstructured]=LimGridDistribX(unstructured,xMax,xMin,maxExcess);

        otherwise
            warning([gridDistrib,' is not a recognised grid distribution Type'])
    end
    
end


function [unstructured]=LimGridDistribX(unstructured,xMax,xMin,maxExcess)
    
    coord=unstructured.vertex.coord;
    x=coord(:,1);
    
    posiExcess=max(x-xMax);
    negiExcess=max(-x+xMin);
    
    %     Dx=xMax-x;
    %     Dx=min(abs(Dx(Dx>1e-10)));
    %     xNorm=(coord(:,1)-xMin)/(xMax-xMin);
    newX=x;
    
    %     DminNx=newX(xNorm<1 & xNorm>0)-1;
    %     DminNx=min(abs(DminNx(DminNx~=0)));
    newX((x>xMax))=(x(x>xMax)-xMax)*maxExcess/posiExcess+xMax;
    newX((x<xMin))=(x(x<xMin)-xMin)*maxExcess/negiExcess+xMin;
    
    coord(:,1)=newX;
    unstructured.vertex.coord=coord;
    
end

function [unstructured]=CosGridDistribX(unstructured,xMax,xMin)
    
    coord=unstructured.vertex.coord;
    x=coord(:,1);
    
    Dx=xMax-x;
    
    Dx=min(abs(Dx(Dx>1e-10)));
    
    xNorm=(coord(:,1)-xMin)/(xMax-xMin);
    newX=(1-cos(xNorm*pi))/2*(xMax-xMin)+xMin;
    DminNx=newX(xNorm<1 & xNorm>0)-1;
    DminNx=min(abs(DminNx(DminNx~=0)));
    
    newX((x>xMax))=(x(x>xMax)-xMax)/Dx*DminNx+xMax;
    
    newX((x<xMin))=(x(x<xMin)-xMin)/Dx*DminNx+xMin;
    
    coord(:,1)=newX;
    unstructured.vertex.coord=coord;
    
end
function [unstructured]=HalfCosGridDistribX(unstructured,xMax,xMin)
    
    coord=unstructured.vertex.coord;
    x=coord(:,1);
    
    Dx=xMax-x;
    
    Dx=min(abs(Dx(Dx>1e-10)));
    
    xNorm=(coord(:,1)-xMin)/(xMax-xMin);
    newX=(1-cos(xNorm*pi))/2*(xMax-xMin)+xMin;
    DminNx=newX(xNorm<1 & xNorm>0)-1;
    DminNx=min(abs(DminNx(DminNx~=0)));
    
    newX((x>xMax))=(x(x>xMax)-xMax)/Dx*DminNx+xMax;
    
    newX((x<xMin))=(x(x<xMin)-xMin)/Dx*DminNx+xMin;
    newX(x>((xMax+xMin)/2))=x(x>((xMax+xMin)/2));
    coord(:,1)=newX;
    unstructured.vertex.coord=coord;
    
end

function [newnum]=BilinearDistrib(xMax,xMin,centre,num)
    
    currCentre = (xMax+xMin)/2;
    
    newnum = ((num-xMin)/(currCentre-xMin)*(centre-xMin)+xMin).*(num<currCentre)...
        +((num-xMax)/(currCentre-xMax)*(centre-xMax)+xMax).*(num>=currCentre);
end

function [unstructured]=BiLinearGridDistribX(unstructured,xMax,xMin,centre)
    
    coord=unstructured.vertex.coord;
    x=coord(:,1);
    
    newX=BilinearDistrib(xMax,xMin,centre,x);
    
    coord(:,1)=newX;
    unstructured.vertex.coord=coord;
    
end

function [unstructured]=Cos2GridDistribX(unstructured,xMax,xMin)
    
    coord=unstructured.vertex.coord;
    x=coord(:,1);
    
    Dx=xMax-x;
    
    Dx=min(abs(Dx(Dx>1e-10)));
    
    xNorm=(coord(:,1)-xMin)/(xMax-xMin);
    newX=((1-cos(xNorm*pi))/2).^1.4*(xMax-xMin)+xMin;
    DminNx=[newX(xNorm<1 & xNorm>0)-1];
    DminNxBack=min(abs(DminNx(DminNx~=0)));
    DminNx=[newX(xNorm<1 & xNorm>0)];
    DminNxFront=min(abs(DminNx(DminNx~=0)));
    
    newX((x>xMax))=(x(x>xMax)-xMax)/Dx*DminNxBack+xMax;
    
    newX((x<xMin))=(x(x<xMin)-xMin)/Dx*DminNxFront+xMin;
    
    coord(:,1)=newX;
    unstructured.vertex.coord=coord;
    
end

function [unstructured]=CosGridDistribY(unstructured)
    % Centred
    distrib=@(y) (sign(y).*(1-cos(y*pi)).^1.7)/2;
    coord=unstructured.vertex.coord;
    y=coord(:,2);
    y1=unique(y);
    [~,I]=min(y1);y1(I)=[];yMin=min(y1);[~,I]=max(y1);y1(I)=[];yMax=max(y1);
    Dy=yMax-y;
    
    Dy=min(abs(Dy(Dy>1e-10)));
    
    yNorm=(y-yMin)/(yMax-yMin);
    newY=distrib(yNorm-0.5)*(yMax-yMin);
    DminNx=[newY(yNorm<1 & yNorm>0)-yMax];
    DminNxBack=min(abs(DminNx(DminNx~=0)));
    DminNx=[newY(yNorm<1 & yNorm>0)-yMin];
    DminNxFront=min(abs(DminNx(DminNx~=0)));
    
    newY((y>yMax))=(y(y>yMax)-yMax)/Dy*DminNxBack+yMax;
    
    newY((y<yMin))=(y(y<yMin)-yMin)/Dy*DminNxFront+yMin;
    
    coord(:,2)=newY;
    unstructured.vertex.coord=coord;
    
end

% Fill functions

function [unstructReshape] = FillNewGrid(unstructReshape, param)
    
    varExtract={'initialfill'};
    [initialfill]=ExtractVariables(varExtract,param);
    switch initialfill{1}
        case {'none','n/a'}
        case 'load'
            [fill]=MatchVoltoShapeGeneral(unstructReshape,initialfill{2});
        otherwise
            
            warning([initialfill{1},' is not a recognised fill type'])
    end
    
    if exist('constrVal', 'var')
        for ii = 1:numel(unstructReshape.cell)
            unstructReshape.cell(ii).fill = constrVal{2}(ii);
        end
    end
    if exist('fill', 'var')
        for ii = 1:numel(unstructReshape.cell)
            unstructReshape.cell(ii).fill = fill(ii);
        end
        
    end
        
end

%% Initialisation Functions

function [unstructured]=InitialisEdgeGrid(param)
    % Main function for the execution of the Subdivision process
    varExtract={'cellGeometry'};
    [cellGeometry]=ExtractVariables(varExtract,param);
    
    switch cellGeometry
        case 'square'
    [unstructured]=Initialisation_Square(param);
    %edgeTemplate=EdgeBuildTemplate(unstructured);
    %unstructured.edge=CreateEdges(unstructured,edgeTemplate);
        case 'triangle'
            [unstructured]=BuildTriangularGrid(param);
    end
end

function [fill]=InputData(typDat)
    global nDim nGridSteps nPadding
    % Creates geometry information data in a cell centred array
    arraySize=ones(1,nDim)*nGridSteps;
    
    switch typDat
        case 'rand'
            fill=randi(2,arraySize)-1;
            [fill]=AddZeroLayer(fill,nPadding);
            sizeIm=size(fill);
            nGridSteps=sizeIm(1);
        case 'cyli'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'wedg'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'foil'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.jpg'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'foillo'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'foillosubdiv'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.bmp'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uni'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uniLo1'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uniLo2'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'uniLo3'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'square'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'low5shapek'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'k',nPadding);
        case 'low5shape'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'vlofoil'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'w',nPadding);
        case 'doubleBody'
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            fill=ImageProcess(imPath,'k',nPadding);
            
        otherwise
            imPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.png'];
            if strcmp(typDat(end-1:end),'_f')
                
                fill=ImageProcess(imPath,'wf',nPadding);
            else
                fill=ImageProcess(imPath,'wc',nPadding);
            end
    end
    
end

function []=WriteCellGridDat(arraySize,cellRef)
    
    fileName=[cd,'\MEX_Function_Directory\MEX_Executables\gridgen\cellgrid.dat'];
    fileName=MakePathCompliant(fileName);
    fID=fopen(fileName,'w');
    
    fprintf(fID,'1\n');
    for ii=1:length(arraySize(:,1))
        fprintf(fID,'%i %i\n',arraySize(ii,1),arraySize(ii,2));
    end
    if (cellRef(1))~=0
        fprintf(fID,'%i\n',length(cellRef(:,1)));
        
        for ii=1:length(cellRef(:,1))
            fprintf(fID,'%i %i ',cellRef(1),ii);
            for jj=1:length(cellRef(ii,:)-1)
                fprintf(fID,'%i ',cellRef(1+jj));
            end
            fprintf(fID,'\n');
        end
    else
        fprintf(fID,'0\n');
    end
    
end

function [unstructured]=Initialisation_Square(param)
    %
    
    global nGridSteps domainBounds
    varExtract={'typDat'};
    [typDat]=ExtractVariables(varExtract,param);
    
    switch typDat
        case 'optimInit'
            [parametrisation,cellRef]=OptimInit(param,false);
        case 'optimRand'
            [parametrisation,cellRef]=OptimInit(param,true);
        otherwise
            parametrisation.fill=InputData(typDat);
            parametrisation.isactive=ones(size(parametrisation.fill));
            cellRef=[0];
    end
    
    arraySize=nGridSteps;
    
    WriteCellGridDat(arraySize,cellRef);
    [unstructReshape]=GridInit_MEX;
    clear GridInit_MEX
    
    for ii=1:length(unstructReshape.cell)
        unstructReshape.cell(ii).fill=parametrisation.fill(ii);
        unstructReshape.cell(ii).isactive=parametrisation.isactive(ii);
    end
    
    for ii=1:length(unstructReshape.edge)
        unstructReshape.edge(ii).cellindex(unstructReshape.edge(ii).cellindex<0)=0;
    end
    
    domMultiplier=(domainBounds(:,2)-domainBounds(:,1))';
    domAdd=domainBounds(:,1)';
    
    for ii=1:length(unstructReshape.vertex)
        unstructReshape.vertex(ii).coord=unstructReshape.vertex(ii).coord.*domMultiplier+domAdd;
    end
    unstructured=ModifReshape(unstructReshape);
    
end

%% Optimisation Input

function [parametrisation,cellRef]=OptimInit(param,isRand)
    
    if ~exist('isRand','var'),isRand=false;end
    
    global nGridSteps
    
    varExtract={'defaultfill','cellLevels','refineCellLvl','defaultCorner','corneractive'};
    [defaultfill,cellLevels,refineCellLvl,defaultCorner,corneractive]=ExtractVariables(varExtract,param);
    
    
    nGridSteps=cellLevels;
    nGridSteps(1,:)=nGridSteps(1,:)+2;
    
    cellRef=refineCellLvl;
    parametrisation.fill=zeros(nGridSteps);
    parametrisation.isactive=zeros(nGridSteps);
    if ~isRand
        parametrisation.fill(2:nGridSteps(1)-1,2:nGridSteps(2)-1)=defaultfill;
    else
        parametrisation.fill(2:nGridSteps(1)-1,2:nGridSteps(2)-1)=...
            rand(size(parametrisation.fill(2:nGridSteps(1)-1,2:nGridSteps(2)-1)));
    end
    parametrisation.isactive(2:nGridSteps(1)-1,2:nGridSteps(2)-1)=true;
    if islogical(corneractive) || (~ischar(corneractive) && (corneractive==0 || corneractive==1))
        
        
        parametrisation.fill([2,nGridSteps(1)-1],[2,nGridSteps(2)-1])=defaultCorner;
        parametrisation.isactive([2,nGridSteps(1)-1],[2,nGridSteps(2)-1])=corneractive;
        
        parametrisation.fill=parametrisation.fill;
    else
        LEcol=2;
        TEcol=nGridSteps(1)-1;
        vertPosRow=floor(nGridSteps(2)/2);
        vertPosRow=vertPosRow+rem(nGridSteps(2),2):vertPosRow+1;
        switch corneractive
            case 'LE'
                parametrisation.isactive(LEcol,:)=false;
                parametrisation.fill(LEcol,:)=0;
                parametrisation.fill(LEcol,vertPosRow)=defaultCorner;
            case 'TE'
                parametrisation.isactive(TEcol,:)=false;
                parametrisation.fill(TEcol,:)=0;
                parametrisation.fill(TEcol,vertPosRow)=defaultCorner;
            case 'LETE'
                parametrisation.isactive([LEcol,TEcol],:)=false;
                parametrisation.fill([LEcol,TEcol],:)=0;
                parametrisation.fill([LEcol,TEcol],vertPosRow)=defaultCorner;
        end
    end
    
    
    
end

%% Treatment of Input Images
function finishedImage=ImageProcess(imPath,imType,nPad,n)
    % Processes images into a valid input to the program
    % impath indicates the background colour: 'k' is black and 'w' is white
    
    global nGridSteps % number of steps in design domain
    if ~exist('n','var'); n=0; end % n is used to pad with ones instead of zeros
    
    if numel(imType)==1
        preProcImage=PreProcImage(imPath);
    elseif strcmp(imType(2),'c')
        
        preProcImage=PreProcImage(imPath);
    elseif strcmp(imType(2),'f')
        
        [preProcImage]=ProcImageFine(imPath);
    end
    
    preProcImage=ProcessType(imType,preProcImage);
    procImage=GradProcessor(preProcImage);
    finishedImage=ResizeImage(procImage);
    [finishedImage]=AddZeroLayer(finishedImage,nPad,n);
    sizeIm=size(finishedImage);
    
    nGridSteps=sizeIm;
    
end

function [preProcImage]=PreProcImage(imPath)
    % Load Image and reduce it to an averaged double array from 0 to 1
    
    preProcImage=imread(MakePathCompliant(imPath));
    imClass=class(preProcImage);
    numBit=str2num(regexprep(imClass,'uint',''));
    preProcImage=mean(preProcImage,3);
    
    %     lvlMax=max(preProcImage(:));
    %     lvlMin=min(preProcImage(:));
    
    lvlMax=2^numBit-1;
    lvlMin=0;
    
    if lvlMin~=lvlMax
        preProcImage=(preProcImage-lvlMin)/(lvlMax-lvlMin);
    else
        preProcImage=zeros(size(preProcImage));
    end
    
end

function [preProcImage]=ProcImageFine(imPath)
    % Load Image and reduce it to an averaged double array from 0 to 1
    % R is Small
    % Green is inter
    % B is big
    imPath=MakePathCompliant(imPath);
    preProcImage=imread(imPath);
    imClass=class(preProcImage);
    numBit=str2num(regexprep(imClass,'uint',''));
    
    preProcImage=double(preProcImage);

    lvlMax=2^numBit;
    lvlMin=0;
    
    if lvlMin~=(lvlMax-1)

        for ii=1:size(preProcImage,3)
            preProcImage(:,:,ii)=preProcImage(:,:,ii)*lvlMax^(ii-1);
        end
        preProcImage=sum(preProcImage,3)/lvlMax^size(preProcImage,3);
    else
        preProcImage=zeros(size(preProcImage));
    end
    
end

function [imageDat]=ProcessType(imType,imageDat)
    % Processes the background of the image (either white or black)
    
    switch imType(1)
        case 'k'
            imBackground=0;
        case 'w'
            imBackground=1;
    end
    
    imageDat=abs(imageDat-imBackground);
    
end

function [finishedImage]=ResizeImage(procImage)
    % Resizes the image to the right size for the execution of the program
    
    finishedImage=TrimZerosImage(procImage);
    %finishedImage=SquareArray(procImage);
    
    finishedImage=finishedImage';
    finishedImage=finishedImage(:,end:-1:1);
end

function [procImage]=GradProcessor(preProcImage)
    % Process gradients to ones or zeros
    
    procImage=(preProcImage);
end

function [array]=TrimZerosImage(array)
    % this function trims complete lines and columns of zeros from the borders
    % of an array
    
    sizArray=size(array);
    nDim=length(sizArray);
    infoArray{nDim}=[];
    
    for iDim=1:nDim
        dims=1:nDim;
        dims(dims==iDim)=[];
        workingArray=array;
        for ii=dims
            workingArray=sum(workingArray,ii);
        end
        infoArray{iDim}=squeeze(workingArray);
    end
    
    for iDim=1:nDim
        
        indexCell=[];
        indexCell{nDim}=[];
        for ii=1:nDim
            indexCell{ii}=':';
        end
        
        candidates=find(infoArray{iDim}==0);
        if ~isempty(candidates)
            dimRmv=[];
            iRmv=1;
            kkLow=0;
            condition=candidates(kkLow+1)==(kkLow+1);
            while condition
                kkLow=kkLow+1;
                dimRmv(iRmv)=kkLow;
                iRmv=iRmv+1;
                condition=(kkLow+1)<=numel(candidates);
                if condition
                    condition=candidates(kkLow+1)==(kkLow+1);
                end
                
            end
            
            kkHigh=sizArray(iDim)+1;
            kk=length(candidates)+1;
            while (candidates(kk-1)==(kkHigh-1))
                kkHigh=kkHigh-1;
                kk=kk-1;
                dimRmv(iRmv)=kkHigh;
                iRmv=iRmv+1;
                if kk==1
                    break
                end
            end
            indexCell{iDim}=dimRmv;
            array(indexCell{:})=[];
        end
    end
    
end

function [array]=AddZeroLayer(array,nPad,n)
    % Adds n layer of zeros to every side of the array
    
    sizArray=size(array);
    sizNew=sizArray+2*nPad;
    
    nDim=length(sizArray);
    for ii=1:nDim
        sizArray=size(array);
        sizePads=sizArray;
        sizePads(ii)=2*nPad+sizePads(ii);
        newArray=zeros(sizePads)+n;
        for jj=1:nDim
            indexCell{jj}=':';
        end
        indexCell{ii}=nPad+1:sizePads(ii)-nPad;
        newArray(indexCell{:})=array;
        array=newArray;
    end
end

%% Triangular grid option

function [unstructured]=BuildTriangularGrid(param)
    varExtract={'ptsDistrib','cellLevels','defaultfill','passDomBounds'};
    [ptsDistrib,cellLevels,defaultfill,passDomBounds]=ExtractVariables(...
        varExtract,param);
    
    pts=DesignPtsForTriangular(cellLevels,passDomBounds,ptsDistrib);
    
    [unstructReshape]=Pts2DelaunayGrid(pts,defaultfill);
    
    
    [unstructured]=ModifReshape(unstructReshape);
    [~,sortOrd]=sort(unstructured.edge.index);
    for ii=fieldnames(unstructured.edge)
        unstructured.edge.(ii{1})=unstructured.edge.(ii{1})(sortOrd,:);
    end
    
end

function pts=DesignPtsForTriangular(cellLevels,passDomBounds,ptsDistrib)
    
    corners=[0 0;0 1; 1 1;1 0];
    corners=[corners;0.5 0;0 0.5; 0.5 1;1 0.5];
    
    switch ptsDistrib
        case 'lhs'
            nLhs=prod(cellLevels)/2;
            pts=lhsdesign(nLhs,2,'criterion','maximin',...
                'iterations',1000);
            [dist]=DistanceFunction(pts,corners);
            pts(find(any(dist<(1/sqrt(nLhs)/4),2)),:)=[];
            
            pts=[pts;corners];
            
        case 'FF'
            [x,y]=meshgrid(0:cellLevels(1),0:cellLevels(2));
            pts=[x(:)/cellLevels(1),y(:)/cellLevels(2)];
        case 'FFStaggered'
            [x,y]=meshgrid(0:cellLevels(1),0:ceil(cellLevels(2)/2));
            pts=[x(:)/cellLevels(1),y(:)/ceil(cellLevels(2)/2)];
            [x,y]=meshgrid(0:cellLevels(1)-1,0:ceil(cellLevels(2)/2)-1);
            pts2=[x(:)/cellLevels(1)+0.5/cellLevels(1),...
                y(:)/ceil(cellLevels(2)/2)+0.5/ceil(cellLevels(2)/2)];
            pts=[pts;pts2];
        case 'rand'
            pts=rand([prod(cellLevels),2]);
            pts=[pts;corners];
        case 'lhsrep'
            isFlip=cellLevels(1)>cellLevels(2);
            nRep=ceil(max(cellLevels)/min(cellLevels));
            nLhs=ceil(prod(cellLevels)/nRep);
            pts=zeros([0 2]);
            cornTemp=zeros([0 2]);
            for ii=1:nRep
                cornTemp(:,2)=cornTemp(:,2)+1;
                pts(:,2)=pts(:,2)+1;
                ptsTemp=lhsdesign(nLhs,2,'criterion','maximin',...
                    'iterations',1000);
                pts=[pts;ptsTemp];
                cornTemp=[cornTemp;corners];
            end
            [dist]=DistanceFunction(pts,cornTemp);
            pts(find(any(dist<(1/sqrt(nLhs)/4),2)),:)=[];
            pts=[pts;cornTemp];
            pts(:,2)=pts(:,2)/nRep;
            if isFlip
                pts=[pts(:,2),pts(:,1)];
            end
            pts=RemoveIdenticalVectors(pts);
        otherwise
            error('unrecognised point distribution')
    end
    
    for ii=1:2
        pts(:,ii)=pts(:,ii)*(passDomBounds(ii,2)-passDomBounds(ii,1))...
            +passDomBounds(ii,1);
    end
    
end

function [dist]=DistanceFunction(pts,ptsref)
    
    dist=zeros(size(pts,1),size(ptsref,1));
    np=size(pts,1);
    for ii=1:size(ptsref,1)
        dist(:,ii)=sqrt(sum((pts-repmat(ptsref(ii,:),[np,1])).^2,2));
    end
    
end

function [gridtri]=Pts2DelaunayGrid(pts,defaultfill)
    
    dtri=delaunayTriangulation(pts);
    [gridtri]=BasicGrid(dtri,pts);
    [gridtri.cell(:).fill]=deal(defaultfill);
    [gridtri]=ZeroBorderCell(gridtri);
end

function [baseGrid]=ZeroBorderCell(baseGrid)
    
    cellCentredGrid=CellCentredGrid(baseGrid);
    
    for ii=1:numel(baseGrid.cell)
        baseGrid.cell(ii).fill=baseGrid.cell(ii).fill*(all(...
            [cellCentredGrid(ii).edge.cellindex]>0));
    end
    
    
end

function [gridtri]=BasicGrid(dtri,pts)
    
    kke=1;
    ndtri2=size(dtri,2);
    gridtri.cell=repmat(struct('index',0,'fill',[],'isactive',true),[size(dtri,1),1]);
    gridtri.vertex=repmat(struct('index',0,'coord',[0 0]),[size(pts,1),1]);
    gridtri.edge=repmat(struct('index',0,'cellindex',[],'vertexindex',[],'orientation',0.5),[size(dtri,1)*3,1]);
    for ii=1:size(pts,1)
        gridtri.vertex(ii).index=ii;
        gridtri.vertex(ii).coord=pts(ii,:);
    end
    for ii=1:size(dtri,1)
        gridtri.cell(ii).index=ii;
        for jj=1:ndtri2
            gridtri.edge(kke).vertexindex=sort([dtri(ii,jj),dtri(ii,mod(jj,ndtri2)+1)]);
            gridtri.edge(kke).cellindex=ii;
            
            kke=kke+1;
        end
    end
    
    
    cellSameEdge=FindIdenticalVector(vertcat(gridtri.edge(:).vertexindex));
    edgeRm=[];
    for ii=1:numel(cellSameEdge)
        gridtri.edge(cellSameEdge{ii}(1)).cellindex=...
            [gridtri.edge(cellSameEdge{ii}).cellindex];
        if numel(gridtri.edge(cellSameEdge{ii}(1)).cellindex)==1
           gridtri.edge(cellSameEdge{ii}(1)).cellindex(2)=0; 
        end
        gridtri.edge(cellSameEdge{ii}(1)).index=ii;
        edgeRm=[edgeRm,cellSameEdge{ii}(2:end)];
    end
    gridtri.edge(edgeRm)=[];
    
end