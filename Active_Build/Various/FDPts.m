function [points]=FDPts(points,n,m)
    
%     dd1=@(pm1,pp1) (pp1(:,2)-pm1(:,2))./(pp1(:,1)-pm1(:,1));
%     dx1c=@(p0,pm1,pp1) ((pp1(:,2)-p0(:,2)).*(p0(:,1)-pm1(:,1))+...
%         (p0(:,2)-pm1(:,2)).*((pp1(:,1)-p0(:,1))))./...
%         ((pp1(:,1)-p0(:,1)).*(p0(:,1)-pm1(:,1)));

    % first derivative with variable stenci size
    dd1c=@(p0,pm1,pp1) ((pp1(:,2)-p0(:,2)).*(p0(:,1)-pm1(:,1)).^2+...
        (p0(:,2)-pm1(:,2)).*((pp1(:,1)-p0(:,1))).^2)./...
        ((pp1(:,1)-p0(:,1)).*(p0(:,1)-pm1(:,1)).^2+...
        (pp1(:,1)-p0(:,1)).^2.*(p0(:,1)-pm1(:,1)));
    
    dd2c=@(p0,pm1,pp1) 2*((pp1(:,2)-p0(:,2)).*(p0(:,1)-pm1(:,1))-...
        (p0(:,2)-pm1(:,2)).*((pp1(:,1)-p0(:,1))))./...
        (abs(pp1(:,1)-p0(:,1)).*(p0(:,1)-pm1(:,1)).^2+...
        (pp1(:,1)-p0(:,1)).^2.*abs(p0(:,1)-pm1(:,1)));
    
    
    integr1=@(p0,pm1) cumsum((p0(:,2)+pm1(:,2))/2.*(p0(:,1)-pm1(:,1)));
    
    if n>0
        ii=1;
        points(:,ii+2)=dd1c(points(:,[1,ii+1]),points([1,1:end-1],[1,ii+1]),...
            points([2:end,end],[1,ii+1]));
    end
    for ii=2:n
        points(:,ii+2)=dd2c(points(:,[1,ii]),points([1,1:end-1],[1,ii]),...
            points([2:end,end],[1,ii]));
    end
    
    for ii=1:m
        % n*sign(ii-1)==0 for ii=1 which allows to use the point y
        % coordinate and not the derivative.
        points(:,n+2+ii)=integr1(points(:,[1,ii+1+n*sign(ii-1)]),...
            points([1,1:end-1],[1,ii+1+n*sign(ii-1)]));
    end
end