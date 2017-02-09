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
    
    
    [analysisCoord,upperLower]=PrepareLoopCoord(loop);
    
    switch aeroClass
        case 'NACA'
            [targCoord]=GenerateNacaCoord(analysisCoord(:,1),upperLower,aeroName);
        case 'lib'
            
            [targCoord]=GenerateLibCoord(analysisCoord(:,1),upperLower,aeroName);
        otherwise
            error('not coded yet') 
    end
    
    switch profileComp
        case 'distance'
            [errorMeasure,modifiedDistance]=CompareProfilesDistance(analysisCoord,targCoord);
        case 'area'
            
            [errorMeasure,modifiedDistance]=CompareProfilesArea(analysisCoord,targCoord);
            
        otherwise
            error('not coded yet') 
    end
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
    
end




%%

function [analysisCoord,upperLower]=PrepareLoopCoord(loop)
    % prepares the loop into the right coordinates
    
    nLoop=numel(loop);
    analysisCoord=vertcat(loop(:).subdivspline);
    
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
        analysisCoord(rmRow,:)=[];
        upperLower(rmRow)=[];
        
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

function [errorMeasure]=CompareProfilesArea(profileCoord,targCoord)
    
    [x0,y0,iout,jout] = intersections(profileCoord(:,1),profileCoord(:,2),...
        targCoord(:,1),targCoord(:,2));
    targArea=abs(CalculatePolyArea(targCoord));
    nX0=numel(x0);
    nP1=size(profileCoord,1);
    nP2=size(targCoord,1);
    if nX0>0
        areaErr=zeros(size(x0));
        areaLength=zeros(size(x0));
        for ii=1:nX0
            % Need to make sure both profiles go in the same direction
            iip1=mod(ii,nX0)+1;
            if ceil(iout(ii))>floor(iout(iip1))
                ind1=[ceil(iout(ii)):nP1,1:floor(iout(iip1))];
            else
                ind1=ceil(iout(ii)):floor(iout(iip1));
            end
            if ceil(jout(ii))>floor(jout(iip1))
                ind2=flip([ceil(jout(ii)):nP2,1:floor(jout(iip1))]);
            else
                ind2=flip(ceil(jout(ii)):floor(jout(iip1)));
            end
            
            actPts=[[x0(ii),y0(ii)]; 
                profileCoord(ind1,:);
                [x0(iip1),y0(iip1)]; 
                targCoord(ind2,:)];
            
            areaErr(ii)=abs(CalculatePolyArea(actPts));
            areaLength(ii)=x0(ii)-x0(iip1);
        end
        errorMeasure=sum(areaErr)/targArea;
        distribError=areaErr./areaLength;
        xPts=[[profileCoord(1,1),reshape(x0,[1,numel(x0)])];
            reshape(x0,[1,numel(x0)]),profileCoord(end,1)]]
        xPts=xPts(:);
        errPts=[distribError(end),distribError;distribError(end),distribError];
        errPts=errPts(:);
        errorDist=[xPts',errPts'];
    else
        
        [errorMeasure]=abs(CalculatePolyArea(profileCoord)-targArea)/targArea;
    end
    
    %error('Compare Profile through area integration has not been coded yet')
end








