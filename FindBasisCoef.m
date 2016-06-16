
function []=FindBasisCoef()
    
    f= @(x,w,c) cos(w*abs(x)+atan(-c/w)).*exp(-c*abs(x));
    f_p= @(x,w,c,p) cos(w*abs(x)+p).*exp(-c*abs(x));
    df= @(x,w,c) (c*cos(w*abs(x)+atan(-c/w))+w*sin(w*abs(x)+atan(-c/w))).*exp(-c*abs(x));
    
    x=linspace(-20,20,10001);
    
    d=2;
    w=1;
    c=0.3;
    figure;
    plot(x,f(x,w,c),x,F_abs(x,w,c),x,df_abs(x,w,c),x,F_nodecay(x,w,c))
    figure
    
    subplot(1,2,1)
    hold on
    plot([-20 20],[0 0],'--')
    subplot(1,2,2)
    hold on
    plot([-20 20],[0 0],'--')
    for ii=1:4
        subplot(1,2,1)
        plot(x,f(x,w,c)+f(x,w,c/ii))
        subplot(1,2,2)
        plot(x,F_abs(x,w,c)+F_abs(x,w,c/ii))
    end
    
    
    %     wvec=linspace(0,2*pi,
    
end

function Y=F_abs(x,w,c)
    F= @(x,w,c) (w*sin(w*x+atan(-c/w))-c*cos(w*x+atan(-c/w))).*exp(-c*x)/(c^2+w^2);
    Fneg=@(x,w,c) -(w*sin(w*-x+atan(-c/w))-c*cos(w*-x+atan(-c/w))).*exp(c*x)/(c^2+w^2);
    
    isNeg=x<0;
    x_abs=abs(x);
    
    %Y=isNeg.*Fneg(x,w,c)+(~isNeg).*(F(x,w,c)+Fneg(0,w,c)-F(0,w,c));
    Y=isNeg.*Fneg(x,w,c)+(~isNeg).*(F(x,w,c));
    (Fneg(0,w,c)-F(0,w,c))
    
end

function Y=F_nodecay(x,w,c)
    F= @(x,w,c) (w*sin(w*x+atan(-c/w))-c*cos(w*x+atan(-c/w)));
    Fneg=@(x,w,c) -(w*sin(w*-x+atan(-c/w))-c*cos(w*-x+atan(-c/w)));
    
    isNeg=x<0;
    x_abs=abs(x);
    
    Y=isNeg.*Fneg(x,w,c)+(~isNeg).*(F(x,w,c));
    
end

function Y=df_abs(x,w,c)
    
    df= @(x,w,c) -(c*cos(w*(x)+atan(-c/w))+w*sin(w*(x)+atan(-c/w))).*exp(-c*(x));
    dfNeg= @(x,w,c) (c*cos(-w*(x)+atan(-c/w))+w*sin(-w*(x)+atan(-c/w))).*exp(c*(x));
    
    isNeg=x<0;
    
    
    Y=isNeg.*dfNeg(x,w,c)+(~isNeg).*(df(x,w,c));
end
