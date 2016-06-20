
function []=FindBasisCoef2()
    
    f= @(x,w,c) cos(w*abs(x)+atan(-c/w)).*exp(-c*abs(x));
    f_p= @(x,w,c,p) cos(w*abs(x)+p).*exp(-c*abs(x));
    df= @(x,w,c) (c*cos(w*abs(x)+atan(-c/w))+w*sin(w*abs(x)+atan(-c/w))).*exp(-c*abs(x));
    
    xmin=-20;
    xmax=20;
    x=linspace(xmin,xmax,10001);
    
    d=1;
    w=pi/(2*d);
    c=w;
    figure;
    [F_absCalc,fLim]=F_abs(x,w,c);
    plot(x,f(x,w,c),x,F_absCalc,x,df_abs(x,w,c))
    box=axis;
    hold on
    plot([xmin 0 0 xmax],[0 0 fLim fLim ],'--')
    plot([xmin xmax],[0 0 ],'--')
    
    for ii=[-d:-2*d:xmin,d:2*d:xmax]
        plot([ii,ii],[-2 2],'--','color',[0.3 0.3 0.3])
    end
    axis(box)
    
    figure
    hold on
    fmat=zeros(size(x));
    kk=1;
    for ii=[-d:-2*d:xmin,d:2*d:xmax]
        fmat(kk,:)=f(x+ii,w,c);
        plot(x,fmat(kk,:));
        kk=kk+1;
    end
    plot(x,(sum(fmat)),'r--')
    %%
    d=pi/(2*2);
    figure
    subplot(1,2,1)
    hold on
    subplot(1,2,2)
    hold on
    %d=1;
    c_pairs=[2,2; 3.211,0.7023;3.452,0.1605;3.251 0.6622;3.291 0.5619;3.271 0.602];
    for ii=1:size(c_pairs,1)
        subplot(1,2,2)
        [Y,lim]=F_abs_double(x,2,c_pairs(ii,1),c_pairs(ii,2));
        l_r=plot(x,Y);
        l_f=plot([xmin 0 0 xmax],[0 0 lim lim ],'--');
        l_f.Color=l_r.Color;
        subplot(1,2,1)
        plot(x,f_double(x,2,c_pairs(ii,1),c_pairs(ii,2)));
    end
    subplot(1,2,2)
    box=axis;
    
    for ii=[-d:-2*d:xmin,d:2*d:xmax]
        plot([ii,ii],[-2 2],'--','color',[0.3 0.3 0.3])
    end
    axis(box)
    figure
    hold on
    fmat=zeros(size(x));
    
    for jj=1:size(c_pairs,1)
        kk=1;
        for ii=[-d:-2*d:xmin,d:2*d:xmax]
            fmat(kk,:)=f_double(x+ii,2,c_pairs(jj,1),c_pairs(jj,2));
            %plot(x,fmat(kk,:));
            kk=kk+1;
        end
        l(jj)=plot(x,(sum(fmat)));
    end
    Decay_Relations(1/pi)
    %     figure
    %
    %     subplot(1,2,1)
    %     hold on
    %     plot([-20 20],[0 0],'--')
    %     subplot(1,2,2)
    %     hold on
    %     plot([-20 20],[0 0],'--')
    %     for ii=1:4
    %         subplot(1,2,1)
    %         plot(x,f(x,w,c)+f(x,w,c/ii))
    %         subplot(1,2,2)
    %         plot(x,F_abs(x,w,c)+F_abs(x,w,c/ii))
    %     end
    
    
    %     wvec=linspace(0,2*pi,
    
end

function []=Decay_Relations(w)
    phaseShift=@(c,w) -c.*sin(atan(-w./c)+pi/2)+w*cos(atan(-w./c)+pi/2);
    figure
    x=linspace(0,2*w,500);
    [X,Y]=meshgrid(x,x);
    Z=phaseShift(X,w)+phaseShift(Y,w);
    C=double(Z>0);
    clc
    h=surf(X,Y,Z,C);
    view(0,90)
    h(1).LineStyle='none';
    
end

function [Y,fLim]=F_abs(x,w,c)
    F= @(x,w,c) (w*sin(w*x+atan(-c/w))-c*cos(w*x+atan(-c/w))).*exp(-c*x)/(c^2+w^2);
    Fneg=@(x,w,c) -(w*sin(w*-x+atan(-c/w))-c*cos(w*-x+atan(-c/w))).*exp(c*x)/(c^2+w^2);
    
    isNeg=x<0;
    x_abs=abs(x);
    
    Y=isNeg.*Fneg(x,w,c)+(~isNeg).*(F(x,w,c)+Fneg(0,w,c)-F(0,w,c));
    %Y=isNeg.*Fneg(x,w,c)+(~isNeg).*(F(x,w,c));
    fLim=(Fneg(0,w,c)-F(0,w,c));
    
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

function Y=f_double(x,w,c1,c2)
    
    f= @(x,w,c1,c2) cos(w*abs(x)+atan(-w/c1)).*exp(-c1*abs(x))...
        +cos(w*abs(x)+atan(-w/c2)).*exp(-c2*abs(x));
    
    Y=f(x,w,c1,c2);
end

function [Y,fLim]=F_abs_double(x,w,c1,c2)
    F= @(x,w,c) (w*sin(w*x+atan(-w/c))-c*cos(w*x+atan(-w/c))).*exp(-c*x)/(c^2+w^2);
    Fneg=@(x,w,c) -(w*sin(w*-x+atan(-w/c))-c*cos(w*-x+atan(-w/c))).*exp(c*x)/(c^2+w^2);
    
    isNeg=x<0;
    x_abs=abs(x);
    
    Y1=isNeg.*Fneg(x,w,c1)+(~isNeg).*(F(x,w,c1)+Fneg(0,w,c1)-F(0,w,c1));
    Y2=isNeg.*Fneg(x,w,c2)+(~isNeg).*(F(x,w,c2)+Fneg(0,w,c2)-F(0,w,c2));
    %Y=isNeg.*Fneg(x,w,c)+(~isNeg).*(F(x,w,c));
    fLim=(Fneg(0,w,c1)-F(0,w,c1))+(Fneg(0,w,c2)-F(0,w,c2));
    Y=Y1+Y2;
end

%% Calculate c2 coeff

function []=BisectCoeff()
    
    
    
    
end