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
        'typDat','loadLogical','isCheckRes','boundstr'};
    [passDomBounds,passGridSteps,passPadding,...
        typDat,loadLogical,isCheckRes,boundstr]=ExtractVariables(varExtract,param);
    
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
        [unstructured]=EdgeProperties(unstructured);
        isEdge=unstructured.edge.(boundstr{1});
        cond=boundstr{3};
        [loop]=OrderSurfaceVertex(unstructured,isEdge,cond);
    else
        datPath=[cd,'\Active_Build','\Sample Geometries\',typDat,'.mat'];
        load(datPath);
    end

    if isCheckRes
        %CheckResults(unstructured,loop)
    end
    disp('Reshaping')
    [unstructReshape]=ModifUnstructured(unstructured);
end


%% Initialisation Functions

function [unstructured]=InitialisEdgeGrid(param)
    % Main function for the execution of the Subdivision process
    
    [unstructured]=Initialisation_Square(param);
    %edgeTemplate=EdgeBuildTemplate(unstructured);
    %unstructured.edge=CreateEdges(unstructured,edgeTemplate);
    
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
            fill=ImageProcess(imPath,'w',nPadding);
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
    
    parametrisation.fill([2,nGridSteps(1)-1],[2,nGridSteps(2)-1])=defaultCorner;
    parametrisation.isactive([2,nGridSteps(1)-1],[2,nGridSteps(2)-1])=corneractive;
    
    parametrisation.fill=parametrisation.fill;
    
    
    
end

%% Treatment of Input Images
function finishedImage=ImageProcess(imPath,imType,nPad,n)
    % Processes images into a valid input to the program
    % impath indicates the background colour: 'k' is black and 'w' is white
    
    global nGridSteps % number of steps in design domain
    if ~exist('n','var'); n=0; end % n is used to pad with ones instead of zeros
    
    preProcImage=PreProcImage(imPath);
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

function [imageDat]=ProcessType(imType,imageDat)
    % Processes the background of the image (either white or black)
    
    switch imType
        case 'k'
            imBackground=0;
        case 'w'
            imBackground=1;
    end
    
    imageDat=abs(imageDat-imBackground);
    
end

function [finishedImage]=ResizeImage(procImage)
    % Resizes the image to the right size for the execution of the program
    
    finishedImage=TrimZeros(procImage);
    %finishedImage=SquareArray(procImage);
    
    finishedImage=finishedImage';
    finishedImage=finishedImage(:,end:-1:1);
end

function [procImage]=GradProcessor(preProcImage)
    % Process gradients to ones or zeros
    
    procImage=(preProcImage);
end

function [array]=TrimZeros(array)
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
