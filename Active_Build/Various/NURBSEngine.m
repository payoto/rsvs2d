function []=NURBSEngine()
    
    knots=[0 0 0 0.45 0.5 0.95 1 1 1];
    pts=[0 0; 0.2 0; 0.5 0.5; 0.8 0;1 0];
    w=[1 1 -1 1 1]';
    deg=3;
    plotPoints= @(points,form) plot(points([1:end],1),points([1:end],2),form);
    u=linspace(0,1,101);
    
    hold on;
    %plotPoints(pts,'s-');
    
    [C]=PlotNURBS(u,pts,knots,w,deg);
    
    plotPoints(C,'-')
end

function [C]=PlotNURBS(u,pts,knots,w,deg)
    
    Z=zeros(size(pts,1),numel(u));
    for ii=1:size(pts,1)
        Z(ii,:) = SplineBasis(ii-1,deg,u,knots);
    end
    
    % One shot calculation of NURBS
    C=(Z'*diag(w)*pts)./repmat((Z'*w),[1,size(pts,2)]);
    
end

function Z = SplineBasis(i,p,u,U)  % from Dom
    % recursive definition of the polynomial bases used for B-Splines and
    % NURBS 
    % (c) Dominic Masters - 2017
    
    i = i+1;
    if p==0
        % Z=zeros(1,length(u));
        for j=1:length(u);
            if j==length(u)
                if u(j)>=U(i) && u(j)<=U(i+1)
                    Z(j)=1;
                else
                    Z(j) = 0;
                end
            else
                if u(j)>=U(i) && u(j)<U(i+1)
                    Z(j)=1;
                else
                    Z(j)=0;
                end
            end
        end
    else
        %         Z1=((u-U(i))/(U(i+p)-U(i))).*N(i-1,p-1,u,U);
        %         Z2=((U(i+p+1)-u)/(U(i+p+1)-U(i+1))).*N((i+1)-1,p-1,u,U);
        %         if isnan(Z1)==1;
        %             if isnan(Z2)==1;
        %                 Z=zeros(size(u));
        %             else
        %                 Z=Z2;
        %             end
        %         else
        %             if isnan(Z2)==1;
        %                 Z=Z1;
        %             else
        %                 Z=Z1+Z2;
        %             end
        %         end
        
        t1n = (u-U(i));
        t1d = (U(i+p)-U(i));
        t2n = (U(i+p+1)-u);
        t2d = (U(i+p+1)-U(i+1));
        
        b1 = SplineBasis(i-1,p-1,u,U);
        b2 = SplineBasis((i+1)-1,p-1,u,U);
        
        if t1d ~= 0
            t1 = t1n/t1d;
        else
            t1 = 0;
        end
        
        if t2d ~= 0
            t2 = t2n/t2d;
        else
            t2 = 0;
        end
        
        Z = t1 .* b1 + t2 .* b2;
        
        
        
    end
end