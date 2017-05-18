function []=RealValuedCalcVariations(a,l)
    
    if l>=0
        sol1(a,l)
    else
        sol2(a,l)
    end
    
end


function []=sol1(a,l)
    
    c1=0;
    lfunc=@(c2,a,constr) (constr-(2*a*c2 + ((2*asin(a*(1/(a^2 + c2^2))^(1/2))...
        + 2*a*(1/(a^2 + c2^2))^(1/2)*(1 - a^2/(a^2 + c2^2))^(1/2))*(a^2 + c2^2))/2));
    kk=0;
    b=2;
    while sign(lfunc(-b,a,l))==sign(lfunc(b,a,l))
        b=b*2;
        kk=kk+1;
    end
    kk
    [c2]=GoldenSection_func0(-b,b,1e-10,{a l},lfunc);
    m= -(1/(a^2 + c2^2))^(1/2);
    
    +sqrt(1 -(m*a-c1)^2)/m+c2
    +sqrt(1 -(m*a+c1)^2)/m+c2
    x=linspace(-a,a,101);
    y=(c2 - (1 - (c1 - m*x).^2).^(1/2)/m);
    yp=-(c1 - m*x)./(1 - (c1 - m*x).^2).^(1/2);
    ypp= m./(1 - (c1 - m*x).^2).^(1/2) - (m*(c1 - m*x).^2)./(1 - (c1 - m*x).^2).^(3/2);
    figure,
    subplot(1,3,1)
    hold on
    plot(x,y)
    subplot(1,3,2)
    hold on
    plot(x,yp)
    subplot(1,3,3)
    hold on
    plot(x,ypp-mean(ypp))
    
end

function []=sol2(a,l)
    
    c1=0;
    lfunc=@(c2,a,constr) (constr-(2*a*c2 - ((2*asin(a*(1/(a^2 + c2^2))^(1/2))...
        + 2*a*(1/(a^2 + c2^2))^(1/2)*(1 - a^2/(a^2 + c2^2))^(1/2))*(a^2 + c2^2))/2));
    b=2;
    kk=0;
    while sign(lfunc(-b,a,l))==sign(lfunc(b,a,l))
        b=b*2;
        kk=kk+1;
    end
    kk
    [c2]=GoldenSection_func0(-b,b,1e-10,{a l},lfunc);
    
    m= (1/(a^2 + c2^2))^(1/2);
   -sqrt(1 -(m*a-c1)^2)/m+c2
  -sqrt(1 -(m*a+c1)^2)/m+c2
    x=linspace(-a,a,101);
    
    
    y=c2 - (1 - (c1 - m*x).^2).^(1/2)/m;
    yp=-(c1 - m*x)./(1 - (c1 - m*x).^2).^(1/2);
    ypp= m./(1 - (c1 - m*x).^2).^(1/2) - (m*(c1 - m*x).^2)./(1 - (c1 - m*x).^2).^(3/2);
    figure,
    subplot(1,3,1)
    hold on
    plot(x,y)
    subplot(1,3,2)
    hold on
    plot(x,yp)
    subplot(1,3,3)
    hold on
    plot(x,ypp-mean(ypp))
    
end


function [np]=GoldenSection_func0(lb,hb,tol,extraIn,func)
    
    gr=(sqrt(5)-1)/2;
    
    c=hb-gr * (hb-lb);
    d=lb+gr * (hb-lb);
    
    while abs(c-d)>tol
        
        fc=abs(func(c,extraIn{:}));
        fd=abs(func(d,extraIn{:}));
        
        if fc<fd
            hb=d;
            
        else
            lb=c;
        end
        
        c=hb-gr * (hb-lb);
        d=lb+gr * (hb-lb);
    end
    
    np=(hb+lb)/2;
    
end