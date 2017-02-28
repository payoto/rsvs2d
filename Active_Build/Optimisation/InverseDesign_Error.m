%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%       Check inverse design match
%        for parametric snakes
%             Alexandre Payot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [errorMeasure,h]=InverseDesign_Error(paramoptim,loop)
    
    varExtract={'aeroClass','aeroName','profileComp'};
    [aeroClass,aeroName,profileComp]=ExtractVariables(varExtract,paramoptim);
    varExtract={'typeLoop'};
    [typeLoop]=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    [analysisCoord,upperLower,targPrep]=PrepareLoopCoord(loop,profileComp,typeLoop);
    
    switch aeroClass
        case 'NACA'
            [targCoord]=GenerateNacaCoord(targPrep(:,1),targPrep(:,2),aeroName);
        case 'lib'
            
            [targCoord]=GenerateLibCoord(targPrep(:,1),targPrep(:,2),aeroName);
        otherwise
            error('not coded yet')
    end
    
    switch profileComp
        case 'distance'
            [errorMeasure,modifiedDistance]=CompareProfilesDistance(analysisCoord,targCoord);
            plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
            h=figure;
            subplot(2,1,1)
            plotPoints(analysisCoord)
            hold on
            plotPoints(targCoord)
            legend('snake points','Target points')
            ax=subplot(2,1,2);
            plot(analysisCoord(:,1),modifiedDistance)
            ax.YScale='log';
        case 'area'
            
            [errorMeasure,modifiedDistance]=CompareProfilesArea(analysisCoord,targCoord);
            plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
            h=figure;
            subplot(2,1,1)
            plotPoints(analysisCoord)
            hold on
            plotPoints(targCoord)
            legend('snake points','Target points')
            ax=subplot(2,1,2);
            plotPoints(modifiedDistance)
            ax.YScale='log';
        otherwise
            error('not coded yet')
    end
    
    
end


%%

function [analysisCoord,upperLower,targPrep]=PrepareLoopCoord(loop,profileComp,typeLoop)
    % prepares the loop into the right coordinates
    
    nLoop=numel(loop);
    analysisCoord=vertcat(loop(:).(typeLoop));
    
    if true % nLoop==1;
        [isCCW]=CCWLoop(analysisCoord);
        if ~isCCW
            analysisCoord=flip(analysisCoord,1);
        end
        [~,iTE]=max(analysisCoord(:,1));
        
        analysisCoord=analysisCoord([iTE:end,1:iTE],:); % repeats the trailing edge at start and end.
        %analysisCoord(:,2)=analysisCoord(:,2)-analysisCoord(1,2); % 0 the y at the trailing edge.
        % TE is start and end of coord list, need LE
        %[dLE,iLE]=max(sum((analysisCoord-ones([size(analysisCoord,1),1])*analysisCoord(1,:)).^2,2));
        [~,iLE]=min(analysisCoord(:,1));
        analysisCoord=analysisCoord([1:iLE,iLE:end],:);
        upperLower=ones(size(analysisCoord(:,1)));
        upperLower(iLE+1:end)=-1;
        %analysisCoord(:,1)=analysisCoord(:,1)-min(analysisCoord(:,1));
        rmRow=find(analysisCoord(:,1)>1 | analysisCoord(:,1)<0);
        
        switch profileComp
            case 'distance'
              analysisCoord(rmRow,:)=[];  
              upperLower(rmRow)=[];
              targPrep=[analysisCoord(:,1),upperLower];
            case 'area'
                x=[linspace(1,0,sum(upperLower>0)),...
                    linspace(0,1,sum(upperLower<0))];
              targPrep=[x',upperLower];
                
        end
        
        
    else
        warning('More than one loop to compare')
    end
    
end



%% Profile extraction

function [nacaCoord]=GenerateLibCoord(x,uplow,airfoilstr)
    % Generates y coordinate for a NACA airfoil.GenerateLibCoord
    % X are the desired x coordinates, uplow is a vector the same size
    
    % xTop and xBot need to be normalised
    
    teps=5.48e-04/2/0.8; % true @ corner=1e-5
    
    
    [airfoilDat,airfoil]=ReadAirfoilData(airfoilstr,'');
    
    xMax=airfoil.func.xMax; %max(xPos);
    xMin=airfoil.func.xMin; %min(xPos);
    x(x>xMax)=xMax;
    x(x<xMin)=xMin;
    
    if any(abs(uplow)~=1)
        error('Vector uplow indicating uppper or lower surface is not well formed')
    end
    y=zeros(size(x));
    y(uplow==1)=airfoil.func.upper(x(uplow==1))+(x(uplow==1)-xMin)/(xMax-xMin)*teps;
    y(uplow==-1)=airfoil.func.lower(x(uplow==-1))-(x(uplow==-1)-xMin)/(xMax-xMin)*teps;
    
    
    if size(x,1)==1
        nacaCoord=[x;y]';
    else
        nacaCoord=[x,y];
    end
    
    
end

function [nacaCoord]=GenerateNacaCoord(x,uplow,nacaStr)
    % Generates y coordinate for a NACA airfoil.
    % X are the desired x coordinates, uplow is a vector the same size
    
    % xTop and xBot need to be normalised
    a4_open=0.1015;
    a4_closed=0.1036;
    
    naca45t=@(x,t,c,xMin,a4,teps)  5*t*c*(0.2969*sqrt((x-xMin)/c)-0.1260*((x-xMin)/c)...
        -0.3516*((x-xMin)/c).^2+0.2843*((x-xMin)/c).^3-a4*((x-xMin)/c).^4)+((x-xMin)/c)*teps;
    
    naca4c=@(x,m,p,c,xMin) [m/p^2*(2*p*((x((x-xMin)<(p*c))-xMin)/c)-((x((x-xMin)<(p*c))-xMin)/c).^2),...
        m/(1-p)^2*((1-2*p)+2*p*((x((x-xMin)>=(p*c))-xMin)/c)-((x((x-xMin)>=(p*c))-xMin)/c).^2)];
    
    teps=5.48e-04/2/0.8; % true @ corner=1e-5
    x(x>1)=1;
    x(x<0)=0;
    
    [ctc,pct,tmax,refFlag]=ReadNacaString(nacaStr);
    
    if numel(nacaStr)==4
        
        
        tDist=naca45t(x',tmax,1,0,a4_closed,teps)';
        [x2,sortOrd]=sort(x);
        [~,sortOrd2]=sort(sortOrd);
        cDist=naca4c(x2',ctc,pct,1,0)';
        cDist=cDist(sortOrd2);
        if any(abs(uplow)~=1)
            warning('Vector uplow indicating uppper or lower surface is not well formed')
        end
        y=cDist+uplow.*tDist;
        if size(x,1)==1
            nacaCoord=[x;y]';
        else
            nacaCoord=[x,y];
        end
    elseif numel(nacaStr)==5
        
        error('Five digits not implemented - need for tabulated data')
    else
        error('Unrecognised naca String length')
        
    end
    
end


%% Error Matching

function [errorMeasure,modifiedDistance]=CompareProfilesDistance(profileCoord,targCoord)
    % Compares point to point distance
    
    multipliers=ones(size(targCoord(:,1)));
    multipliers(targCoord(:,1)<0.2*max(targCoord(:,1)))=2;
    
    modifiedDistance=sqrt(sum((profileCoord-targCoord).^2,2)).*multipliers;
    
    errorMeasure.sum=sum(modifiedDistance);
    errorMeasure.mean=mean(modifiedDistance);
    errorMeasure.std=std(modifiedDistance);
    errorMeasure.max=max(modifiedDistance);
    errorMeasure.min=min(modifiedDistance);
    
end

function [errorMeasure,areaDistrib]=CompareProfilesArea(profileCoord,targCoord)
    
    [x0,y0,iout,jout] = intersections(profileCoord([1:end,1],1),profileCoord([1:end,1],2),...
        targCoord([1:end,1],1),targCoord([1:end,1],2));
    rmv=isnan(iout) | isnan(jout);
    iout(rmv)=[];
    x0(rmv)=[];
    y0(rmv)=[];
    jout(rmv)=[];
    
    [iout,sortIndex]=sort(iout);
    x0=x0(sortIndex);
    y0=y0(sortIndex);
    jout=jout(sortIndex);
    targArea=abs(CalculatePolyArea(targCoord));
    nX0=numel(x0);
    nP1=size(profileCoord,1);
    nP2=size(targCoord,1);
    
    if nX0>0
        areaErr=zeros(size(x0))';
        areaLength=zeros(size(x0))';
        areaPosx=zeros(size(x0))';
        areaPosXmin=zeros(size(x0))';
        areaPosXmax=zeros(size(x0))';
        for ii=1:nX0
            % Need to make sure both profiles go in the same direction
            iip1=mod(ii,nX0)+1;
            if floor(iout(ii))~=floor(iout(iip1))
                if ceil(iout(ii))>floor(iout(iip1))
                    ind1=[ceil(iout(ii)):nP1,1:floor(iout(iip1))];
                else
                    ind1=ceil(iout(ii)):floor(iout(iip1));
                end
            else
                ind1=[];
            end
            if floor(jout(ii))~=floor(jout(iip1))
                if ceil(jout(ii))>floor(jout(iip1))
                    ind2=flip([ceil(jout(ii)):nP2,1:floor(jout(iip1))]);
                else
                    ind2=flip(ceil(jout(ii)):floor(jout(iip1)));
                end
            else
                ind2=[];
            end
            
            actPts=[[x0(ii),y0(ii)];
                profileCoord(ind1,:);
                [x0(iip1),y0(iip1)];
                targCoord(ind2,:)];
            areaErr(ii)=abs(CalculatePolyArea(actPts));
            areaLength(ii)=max(actPts(:,1))-min(actPts(:,1));
            areaPosx(ii)=(min(actPts(:,1))+max(actPts(:,1)))/2;
            areaPosXmin(ii)=min(actPts(:,1));
            areaPosXmax(ii)=max(actPts(:,1));
        end
        
        errorMeasure.sum=sum(areaErr)/targArea;
        errorMeasure.mean=mean(areaErr)/targArea;
        errorMeasure.std=std(areaErr)/targArea;
        errorMeasure.max=max(areaErr)/targArea;
        errorMeasure.min=min(areaErr)/targArea;
        
        areaDistrib=[min(areaPosXmin),areaPosx,areaPosXmax;
            0,areaErr./areaLength*2,-areaErr./areaLength*2];
        [areaDistrib(1,:),iSortDistrib]=sort(areaDistrib(1,:));
        eqInd=find(abs(areaDistrib(1,1:end-1)-areaDistrib(1,2:end))<1e-10);
        areaDistrib(2,:)=(areaDistrib(2,iSortDistrib));
        areaDistrib(2,eqInd)=areaDistrib(2,eqInd)+areaDistrib(2,eqInd+1);
        areaDistrib(:,eqInd+1)=[];
        areaDistrib(2,:)=cumsum(areaDistrib(2,:));
        areaDistrib(:,find(areaDistrib(2,:)<1e-18))=[];
        areaDistrib=areaDistrib';
        %         distribError=areaErr./areaLength;
        %         xPts=[[profileCoord(1,1),reshape(x0,[1,numel(x0)])];
        %             [reshape(x0,[1,numel(x0)]),profileCoord(end,1)]];
        %         xPts=xPts(:);
        %         errPts=[distribError(end),distribError;distribError(end),distribError];
        %         errPts=errPts(:);
        %         errorDist=[xPts',errPts'];
    else
        
        [errorMeasure]=(abs(CalculatePolyArea(profileCoord))-targArea)/targArea;
        areaDistrib=[(min(profileCoord(:,1))+min(profileCoord(:,1)))/2,errorMeasure];
    end
    
    %error('Compare Profile through area integration has not been coded yet')
end








