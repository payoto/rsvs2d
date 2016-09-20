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




function [errorMeasure]=InverseDesign_Error(paramoptim,loop)
    
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
            [errorMeasure]=CompareProfilesDistance(analysisCoord,targCoord);
            
        otherwise
            error('not coded yet') 
    end
%     plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
%     figure,
%     plotPoints(analysisCoord)
%     hold on
%     plotPoints(targCoord)
    
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
        analysisCoord(:,2)=analysisCoord(:,2)-analysisCoord(1,2); % 0 the y at the trailing edge.
        % TE is start and end of coord list, need LE
        %[dLE,iLE]=max(sum((analysisCoord-ones([size(analysisCoord,1),1])*analysisCoord(1,:)).^2,2));
        [~,iLE]=min(analysisCoord(:,1));
        analysisCoord=analysisCoord([1:iLE,iLE:end],:);
        upperLower=ones(size(analysisCoord(:,1)));
        upperLower(iLE+1:end)=-1;

        
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
    naca45t=@(x,t,a4)  5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-a4*x.^4);
    naca4c=@(x,m,p) [m*x(x<p)/p^2.*(2*p-x(x<p)) ; m*(1-x(x>=p))/(1-p)^2.*(1+x(x>=p)-2*p)];
    
    
    if numel(nacaStr)==4
        m=str2num(nacaStr(1))/100;
        p=str2num(nacaStr(2))/10;
        tmax=str2num(nacaStr(3:4))/100;
        
        tDist=naca45t(x,tmax,a4_closed);
        cDist=naca4c(x,m,p);
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
        m=str2num(nacaStr(1))/100;
        p=str2num(nacaStr(2))/10;
        reflexFlag=str2num(nacaStr(3));
        tmax=str2num(nacaStr(4:5))/100;
        error('Five digits not implemented - need for tabulated data')
    else
        error('Unrecognised naca String length')
        
    end
    
end

function []=GenerateUUICCoord()
    
    
    
    
end


%% Error Matching

function [errorMeasure]=CompareProfilesDistance(profileCoord,targCoord)
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

