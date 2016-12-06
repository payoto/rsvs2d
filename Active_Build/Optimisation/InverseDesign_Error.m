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
            
        otherwise
            error('not coded yet') 
    end
    
    switch profileComp
        case 'distance'
            [errorMeasure,modifiedDistance]=CompareProfilesDistance(analysisCoord,targCoord);
            
        otherwise
            error('not coded yet') 
    end
    plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
    h=figure;
    subplot(2,1,1)
    plotPoints(analysisCoord)
    hold on
    plotPoints(targCoord)
    legend('snake points','NACA points')
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
    multipliers(targCoord(:,1)<0.2*max(targCoord(:,1)))=1;
    
    modifiedDistance=sqrt(sum((profileCoord-targCoord).^2,2)).*multipliers;
    
    errorMeasure.sum=sum(modifiedDistance);
    errorMeasure.mean=mean(modifiedDistance);
    errorMeasure.std=std(modifiedDistance);
    errorMeasure.max=max(modifiedDistance);
    errorMeasure.min=min(modifiedDistance);
    
end

function [errorMeasure]=CompareProfilesArea(profileCoord,targCoord)
    
    error('Compare Profile through area integration has not been coded yet')
end

