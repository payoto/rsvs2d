
function []=SymbolicCalcVariations(exec)
    
    for ii=exec
        switch ii
            case 1
                SymbolicCalc();
            case 2
                DerivativeOp();
            case 3
                HigherDerivativeOrders();
            case 4
                HigherDerivativeOrders2();
        end
    end
            
        
    
end


%% Generate Solutions
function []=SymbolicCalc()
    syms a c1 c2 m l x
    syms f1 f2 f3 y
    f1=sqrt(1-(m*a-c1)^2)/m+c2;
    f2=sqrt(1-(m*a+c1)^2)/m+c2;
    % f1=1-(m*a-c1)^2-m^2*c2^2;
    % f2=1-(m*a+c1)^2-m^2*c2^2;
    f3=2*a*c2+1/(2*m)*((m*a-c1)*sqrt(1-(c1-a*m)^2)-asin(c1-m*a)+(m*a+c1)*sqrt(1-(c1+a*m)^2)+asin(c1+m*a))-l;
    y=sqrt(1-(m*x-c1)^2)/m+c2;
    assume(a>0);
    solSys2=solve([f1 f2 f3],'returnconditions',true);
    % assume(l,'real');
    assume(c2,'real');
    assume(c1,'real');
    assume(m,'real');
    assume(a,'real');
    %assume(l~=0);
    solSys=solve([f1 f2 f3],'returnconditions',true);
    
    
    %% Output Latex code
    
    [str]=MakeLatexEquationCode(y,'y(x)=','');
    [str]=[str;MakeLatexEquationCode(f1,'\mathrm{end\ condition\ 1:\ }','=0')];
    [str]=[str;MakeLatexEquationCode(f2,'\mathrm{end\ condition\ 2:\ }','=0')];
    [str]=[str;MakeLatexEquationCode(f3,'\mathrm{Integral\ condition\ :\ }','=0')];
    
    
    %sol No assumption
    
    
    %sol assumption
    
    
    
    [str2]=[MakeLatexEquationCode(solSys2.c1,'c1=','')];
    [str2]=[str2;MakeLatexEquationCode(solSys2.l,'l=','')];
    [str2]=[str2;MakeLatexEquationCode(solSys2.m,'m=','')];
    [str2]=[str2;MakeLatexEquationCode(solSys2.conditions,'\mathrm{Condition\ :\ }','')];
    
    %% Test conditions
    
    syms cond1(a,c2) cond2(a,c2)
    
    cond1(a,c2)=solSys2.conditions(1);
    cond2(a,c2)=solSys2.conditions(2);
    
    % cond1R=@(a,c2) (a^3*(2*asin(((4*a^2 + c2^4)^(1/2) - c2^2)/(2*a)) + ((1 - ((4*a^2 + c2^4)^(1/2) - c2^2)^2/(4*a^2))^(1/2)*((4*a^2 + c2^4)^(1/2) - c2^2))/a) + c2*((4*a^2 + c2^4)^(1/2) - c2^2)^2 ~= 0 | (4*a^2 + c2^4)^(1/2) == c2^2) & in(((4*a^2 + c2^4)^(1/2) - c2^2)/a^2, 'real');
    % cond2R=@(a,c2) in(((4*a^2 + c2^4)^(1/2) + c2^2)/a^2, 'real') & ((4*a^2 + c2^4)^(1/2) + c2^2 == 0 | a^3*(2*asin(((4*a^2 + c2^4)^(1/2) + c2^2)/(2*a)) + ((1 - ((4*a^2 + c2^4)^(1/2) + c2^2)^2/(4*a^2))^(1/2)*((4*a^2 + c2^4)^(1/2) + c2^2))/a) ~= c2*((4*a^2 + c2^4)^(1/2) + c2^2)^2);
    
    cond1R=@(a,c2)  0 < a.^2 + c2.^2 & 4.*a.*c2 ~= (2.*asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + 2.*a.*(1./(a.^2 + c2.^2)).^(1./2).*(1 - a.^2./(a.^2 + c2.^2)).^(1./2)).*(a.^2 + c2.^2);
    cond2R=@(a,c2)  4.*a.*c2 + (2.*asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + 2.*a.*(1./(a.^2 + c2.^2)).^(1./2).*(1 - a.^2./(a.^2 + c2.^2)).^(1./2)).*(a.^2 + c2.^2) ~= 0 & 0 < a.^2 + c2.^2;
    
    
    constr1=@(a,c2)  2.*a.*c2 - (a.^2 + c2.^2).*(asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + (a.*abs(c2))./(a.^2 + c2.^2));
    constr2=@(a,c2)  2.*a.*c2 + (a.^2 + c2.^2).*(asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + (a.*abs(c2))./(a.^2 + c2.^2));
    
    
    xa=linspace(0.01,1,101);
    yc2=linspace(-20,5,101);
    
    [Xa,Yc2]=meshgrid(xa,yc2);
    
    iscond1=(cond1R(Xa,Yc2));
    iscond2=(cond2R(Xa,Yc2));
    constrval1=constr1(Xa,Yc2);
    constrval2=constr2(Xa,Yc2);
    
    figure
    subplot(1,2,1)
    surf(Xa,Yc2,real(imag(constrval2)~=0))
    view(0,90)
    subplot(1,2,2)
    s(1)=surf(Xa,Yc2,real(constrval2));
    hold on
    s(2)=surf(Xa,Yc2,real(constrval2)*0);
    s(2).FaceColor=[1 0 0 ];
    s(2).LineStyle='none';
    %%
    
    assume(l,'clear')
    assume(c1,'clear')
    assume(c2,'clear')
    assume(m,'clear')
    assume(a,'clear')
    
    
    %%
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    a=1;
    c1=0;
    c2=[-1 -1+1e-6];
    l=2.*c2 + (2*asin((1./(c2.^2 + 1)).^(1/2)) + 2*(1./(c2.^2 + 1)).^(1/2).*(1 - 1./(c2.^2 + 1)).^(1/2))/(2*(1./(c2.^2 + 1)).^(1/2));
    
    m=(1./(a.^2 + c2.^2)).^(1/2);
    
    x=linspace(-a,a,500);
    
    %y=c2 + (1 - (c1 - m.*x).^2).^(1/2)./m;
    figure,
    ax(3)=subplot(2,4,1);
    plot(real(c2),real(l),real(c2),real(m));
    hold on
    ax(4)=subplot(2,4,5);
    plot(imag(c2),imag(l),imag(c2),imag(m));
    hold on
    ax(1)=subplot(1,4,2);
    hold on
    ax(2)=subplot(1,4,3);
    hold on
    ax(5)=subplot(1,4,4);
    hold on
    for ii=1:numel(c2)
        y=c2(ii) + (1 - (c1 - m(ii)*x).^2).^(1/2)./m(ii);
        plot(ax(1),x,real(y))
        plot(ax(2),x,imag(y))
        plot(ax(5),x,abs(y))
        integData=integr(x,y);
        lReal(ii)=integData(end);
    end
    plot(ax(3),real(c2),real(lReal));
    plot(ax(4),imag(c2),imag(lReal));
end


function []=DerivativeOp()
    %% Comparison of dy/dyc to a parabola
    % This shows that as yc tends to infinity (ie that the profile tends to a
    % straight line the profile tends to a parabola)
    
    dydyc=@(yc,x,a) (sqrt(a^2/yc^2-x.^2/yc^2+1)-1)./(sqrt(a^2/yc^2-x.^2/yc^2+1));
    ypar=@(x,a) -(x+a).*(x-a);
    a=1;
    x=linspace(-a,a,501);
    
    
    h=figure;
    h.Name='Limit of d(y)/d(yc) as yc tends to infinity';
    listYc=[1e7 1e6 1e5 1e4 1e3 1e2 1e1];
    listYc=[0.001 0.01 0.1 1 10 100 1000 10^3.7];
    subplot(1,2,1), hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,dydyc(ii,x,a)/max(abs(dydyc(ii,x,a))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    lineD(kk)=plot(x,(ypar(x,a)),'k--','DisplayName','parabola');
    legend(lineD)
    
    clear lineD
    subplot(1,2,2), hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,(log10(abs(ypar(x,a)-dydyc(ii,x,a)/max(abs(dydyc(ii,x,a)))))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    ii=10^(3.8);
    lineD(kk)=plot(x,(log10(abs(ypar(x,a)-dydyc(ii,x,a)/max(abs(dydyc(ii,x,a)))))),'k');
    lineD(kk).DisplayName=num2str(ii,'%.1e');
    legend(lineD)
    
    %% Calculating the differences in the coordinate reference frame of the profile locally
    
    yfunc=@(x,yc,a) -(yc-sqrt(yc^2-(x.^2-a^2)));
    dydyc=@(x,yc,a) -(yc~=0)*(sqrt(a^2-x.^2+yc.^2)-yc)./(sqrt(a^2-x.^2+yc.^2)+(yc==0))-0.5*(yc==0);
    dydx=@(x,yc,a) -x./(sqrt(yc.^2-(x.^2-a^2)));
    
    % none of these are right
    dydnyc1=@(x,yc,a) -(yc-sqrt(yc^2-x.^2+a^2))./sqrt(6*yc^2+8*a^2-4*x.^2-2*yc*sqrt(yc^2-4*x.^2+4*a^2));
    dydnyc2=@(x,yc,a,dydyc,dydx) -dydyc(x,yc,a)./sqrt(dydyc(x,yc,a).^2+dydx(x,yc,a).^2+1);
    dydnyc3=@(x,yc,a,dydyc,dydx) -dydyc(x,yc,a)./sqrt(dydyc(x,yc,a).^2+dydx(x,yc,a).^2+1).^2;
    % norm of difference between 3D gradient vec and XY gradient (normalised)
    dydnyc4=@(x,yc,a,dydyc,dydx) sqrt(sum(((cat(3,-dydx(x,yc,a),-dydyc(x,yc,a),ones(size(x)))...
        ./repmat(sqrt(sum(cat(3,-dydx(x,yc,a),-dydyc(x,yc,a),ones(size(x))).^2,3)),[1 1 3]))-...
        (cat(3,-dydx(x,yc,a),zeros(size(x)),ones(size(x)))...
        ./repmat(sqrt(sum(cat(3,-dydx(x,yc,a),zeros(size(x)),ones(size(x))).^2,3)),[1 1 3]))).^2,3));
    
    dydnyc5=@(x,yc,a,dydyc,dydx) sqrt(sum((cat(3,-dydx(x,yc,a),-dydyc(x,yc,a),ones(size(x)))...
        -cat(3,-dydx(x,yc,a),zeros(size(x)),ones(size(x)))...
        ).^2,3));
    
    % Could be this
    dydnyc6=@(x,yc,a,dydyc,dydx) -dydyc(x,yc,a)./sqrt(dydx(x,yc,a).^2+1);
    
    normVec=@(pts) ([0 -1;1 0]*pts')';
    ypar=@(x,a) -(x+a).*(x-a);
    a=1;
    %listYc=[0  0.01 0.1 1 10 100 250];
    listYc=[0 0.01 1 100 10000];
    x=(cos(linspace(-a,a,4001)/a*pi)*a)';
    
    kk=1;
    areaDistrib=cell(1,numel(listYc));
    for ii=listYc
        [~,areaDistrib{kk}]=NormalDistance([x,yfunc(x,ii,a)],[x,yfunc(x,ii+1e-8*(1+ii),a)]);
        kk=kk+1;
    end
    
    %%
    h(1)=figure('Name','Finite Difference Normal response');
    subplot(1,2,1)
    hold on
    kk=1;
    
    for ii=listYc
        lineD(kk)=plot(areaDistrib{kk}(:,1),areaDistrib{kk}(:,2)/max(abs(areaDistrib{kk}(:,2))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    lineD(kk)=plot(x,(ypar(x,a)),'k--','DisplayName','1-x^2');kk=kk+1;
    lineD(kk)=plot(x,yfunc(x,0,a),'k-.','DisplayName','y(x,0,1)');
    legend(lineD)
    clear lineD
    subplot(1,4,3), hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,(log10(abs(ypar(x,a)-areaDistrib{kk}(:,2)/max(abs(areaDistrib{kk}(:,2)))))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    clear lineD
    subplot(1,4,4), hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,(log10(abs(yfunc(x,0,a)-areaDistrib{kk}(:,2)/max(abs(areaDistrib{kk}(:,2)))))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    
    
    h(2)=figure('Name','Analytical Normal response');
    subplot(1,2,1), hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,dydnyc6(x,ii,a,dydyc,dydx)/dydnyc6(0,ii,a,dydyc,dydx));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    lineD(kk)=plot(x,(ypar(x,a)),'k--','DisplayName','1-x^2');kk=kk+1;
    lineD(kk)=plot(x,yfunc(x,0,a),'k-.','DisplayName','y(x,0,1)');
    subplot(1,4,3), hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,(log10(abs(ypar(x,a)-dydnyc6(x,ii,a,dydyc,dydx)/dydnyc6(0,ii,a,dydyc,dydx)))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    clear lineD
    subplot(1,4,4), hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,(log10(abs(yfunc(x,0,a)-dydnyc6(x,ii,a,dydyc,dydx)/dydnyc6(0,ii,a,dydyc,dydx)))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    h(3)=figure('Name','Comparison Analytical to Finite');
    
    hold on;
    kk=1;
    for ii=listYc
        lineD(kk)=plot(x,(log10(abs(dydnyc6(x,ii,a,dydyc,dydx)/dydnyc6(0,ii,a,dydyc,dydx)...
            -areaDistrib{kk}(:,2)/max(abs(areaDistrib{kk}(:,2)))))));
        lineD(kk).DisplayName=num2str(ii,'%.1e');
        kk=kk+1;
    end
    %% 3D plots
    
    yfunc=@(x,yc,a) -(yc-sqrt(yc.^2-(x.^2-a^2)));
    dydyc=@(x,yc,a) -(sqrt(a^2./yc.^2-x.^2./yc.^2+1)-1)./(sqrt(a^2./yc.^2-x.^2./yc.^2+1));
    dydx=@(x,yc,a) -x./(sqrt(yc.^2-(x.^2-a^2)));
    yNorm=@(x,yc,a,dydyc,dydx) cat(3,-dydx(x,yc,a),-dydyc(x,yc,a),ones(size(x)));
    yNorm2=@(x,yc,a) cat(3,x,(sqrt(a^2-x.^2+yc.^2)-yc),sqrt(a^2-x.^2+yc.^2));
    a=1;
    x=linspace(-a,a,101);
    yc=logspace(-3,1,50);
    return;
    [X,YC]=meshgrid(x,yc);
    
    YVAL=yfunc(X,YC,a);
    YVAL = YVAL./repmat(max((abs(YVAL)),[],2),[1 size(X,2)]);
    normVector=yNorm(X,YC,a,dydyc,dydx);
    %norVal=yNorm2(X,YC,a);
    unitNormal=normVector./repmat(sqrt(sum(normVector.^2,3)),[1 1 3]);
    normVector=normVector./repmat(max(sqrt(sum(normVector.^2,3)),[],2),[1 size(X,2) 3]);
    h=figure;
    hold on
    s(1)=surf(X,YC,YVAL);
    s(2)=quiver3(X(:,:),YC(:,:),YVAL(:,:),normVector(:,:,1),normVector(:,:,2),normVector(:,:,3));
    %%
    h(2)=figure;
    for ii=1:3
        ax(ii)=subplot(1,3,ii);
        s2(ii)=surf(X(:,:),YC(:,:),normVector(:,:,ii));
        ax(ii).YScale='log';
    end
    h(3)=figure;
    for ii=1:3
        ax(ii)=subplot(1,3,ii);
        s3(ii)=surf(X(:,:),YC(:,:),normVector(:,:,ii)./repmat(max((abs(normVector(:,:,ii))),[],2),[1 size(X,2)]));
        ax(ii).YScale='log';
    end
    h(4)=figure;
    for ii=1:3
        ax(ii)=subplot(1,3,ii);
        s3(ii)=surf(X(:,:),YC(:,:),unitNormal(:,:,ii));
        ax(ii).YScale='log';
    end
    h(6)=figure;
    for ii=1:3
        ax(ii)=subplot(1,3,ii);
        s3(ii)=surf(X(:,:),YC(:,:),unitNormal(:,:,ii)./repmat(max(abs(unitNormal(:,:,ii)),[],2),[1 size(unitNormal,2)]));
        ax(ii).YScale='log';
    end
    %%
    h(7)=figure;
    unitNormalXY=unitNormal;
    unitNormalXY(:,:,2)=0;
    unitNormalXY=unitNormalXY./repmat(sqrt(sum(unitNormalXY.^2,3)),[1 1 3]);
    crossNorm=(unitNormal-unitNormalXY);
    crossNorm2=sqrt(sum(crossNorm.^2,3));
    for ii=1:3
        ax(ii)=subplot(1,3,ii);
        s3(ii)=surf(X(:,:),YC(:,:),crossNorm(:,:,ii));
        ax(ii).YScale='log';
    end
    h(8)=figure;
    for ii=1:3
        ax(ii)=subplot(1,3,ii);
        s3(ii)=surf(X(:,:),YC(:,:),crossNorm(:,:,ii)./repmat(max(abs(crossNorm(:,:,ii)),[],2),[1 size(crossNorm,2)]));
        ax(ii).YScale='log';
    end
    h(9)=figure;
    surf(X(:,:),YC(:,:),crossNorm2./repmat(max(abs(crossNorm2),[],2),[1 size(crossNorm,2)]));
    %%
    h(5)=figure
    dfdycdx=@(x,yc,a) 2*x.*yc./(yc.^2-4*(x.^2-a^2)).^(3/2);
    gradvec=dydyc(X,YC,a).^2+dydx(X,YC,a).^2; %dfdycdx(X,YC,a);
    gradvec=gradvec./repmat(max((abs(gradvec)),[],2),[1 size(X,2)]);
    surf(X,YC,gradvec);
    
end


function []=HigherDerivativeOrders()
    %% Test higher orders
    
    eq{1,1}=@(x,a,c0,c1) sqrt(a.^2-(c0.*a-x).^2)+c1;
    eq{1,2}=@(x,a,c0,c1) -sqrt(a.^2-(c0.*a-x).^2)+c1;
    eq{2,1}=@(x,a,c0,c1,c2) a.^2./(2).*(acos(1./a.*(c0-x))...
        -(1./a.*(c0-x)).*sqrt(1-(1./a.*(c0-x)).^2))+c1.*x+c2;
    eq{2,2}=@(x,a,c0,c1,c2) -a.^2./(2).*(acos(1./a.*(c0-x))...
        -(1./a.*(c0-x)).*sqrt(1-(1./a.*(c0-x)).^2))+c1.*x+c2;
    eq{3,1}=@(x,a,c0,c1,c2,c3)  -a.^3./(6).*(sqrt(1-(1./a.*(c0-x)).^2).^3+sqrt(1-(1./a.*(c0-x)).^2)...
        -acos(1./a.*(c0-x)).*(1./a.*(c0-x)))+1/2*c1.*x.^2+c2.*x+c3;
    eq{3,2}=@(x,a,c0,c1,c2,c3)  -1*-a.^3./(6).*(sqrt(1-(1./a.*(c0-x)).^2).^3+sqrt(1-(1./a.*(c0-x)).^2)...
        -acos(1./a.*(c0-x)).*(1./a.*(c0-x)))+1/2*c1.*x.^2+c2.*x+c3;
    
    c1eq2=@(c0,a,ym,yp) (ym-(yp-pi*a.^2/2-ym./(c0+a))./(1-1./(c0+a)))./(c0-a);
    c2eq2=@(c0,a,ym,yp) (yp-pi*a.^2/2-ym./(c0+a))./(1-1./(c0+a));
    
    c0=[0 0 0];% defines from 0 to 1
    lm=[1 1 1]; % defines the width
    c1=[1 0 0];
    c2=[0 0 0];
    c3=[0 0 0];
    
    
    
    for ii=1:3
        
        x=linspace(-lm(ii)+c0(ii),lm(ii)+c0(ii),501);
        
        switch ii
            case 1
                subplot(1,3,ii),hold on
                plot(x,eq{ii,1}(x,lm(ii),c0(ii),c1(ii)))
            case 2
                subplot(1,3,2),hold on
                plot(x,eq{ii,1}(x,lm(ii),c0(ii),c1(ii),c2(ii)))
            case 3
                subplot(1,3,3),hold on
                plot(x,eq{ii,1}(x,lm(ii),c0(ii),c1(ii),c2(ii),c3(ii)))
        end
    end
    
    %% analysis of 2nd order
    % second order lacks sufficient degrees of freedom to satisfy area
    % constraints cannot have end conditions and constraint
    c1eq2=@(c0,a,ym,yp) (ym-(yp-pi*a.^2/2-ym.*(c0+a)./(c0-a))./(1-(c0+a)./(c0-a)))./(c0-a);
    c2eq2=@(c0,a,ym,yp) (yp-pi*a.^2/2-ym.*(c0+a)./(c0-a))./(1-(c0+a)./(c0-a));
    
    ym=0;
    yp=1;
    
    c0=0;
    a=[0.5 1 2];
    
    c1=c1eq2(c0,a,ym,yp);
    c2=c2eq2(c0,a,ym,yp);
    
    hold on
    for ii=1:numel(a)
        x=linspace(-a(ii)+c0,a(ii)+c0,501);
        plot(x,eq{2,1}(x,a(ii),c0,c1(ii),c2(ii)))
    end
    
    %% analysis of 3rd order
    % Function Shapes
    
    c2eq3{1}=@(c0,a,ym,yp,c2) (yp-ym+c2/2.*(-(c0+a).^2+(c0-a).^2)+a.^3*pi/6)./(2*a);
    c3eq3{1}=@(c0,a,ym,yp,c2) ym-c2/2.*(c0-a).^2-(c0-a).*(yp-ym+c2/2.*(-(c0+a).^2+(c0-a).^2)+a.^3*pi/6)./(2*a);
    
    ym=0;
    yp=1;
    
    c0=0;
    c1=[0 -1 -2 -3];
    
    a=ones(size(c1));
    
    c2=c2eq3{1}(c0,a,ym,yp,c1);
    c3=c3eq3{1}(c0,a,ym,yp,c1);
    figure
    subplot(1,2,1)
    hold on
    clear l
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501);
        l(ii)=plot(x,eq{3,1}(x,a(ii),c0,c1(ii),c2(ii),c3(ii)));
        l(ii).DisplayName=['+ c1:',num2str(c1(ii),'%.1e')];
    end
    
    % analysis of 3rd order (bis)
    %
    
    c2eq3{2}=@(c0,a,ym,yp,c2) (yp-ym+c2/2.*(-(c0+a).^2+(c0-a).^2)-a.^3*pi/6)./(2*a);
    c3eq3{2}=@(c0,a,ym,yp,c2) ym-c2/2.*(c0-a).^2-(c0-a).*(yp-ym+c2/2.*(-(c0+a).^2+(c0-a).^2)-a.^3*pi/6)./(2*a);
    
    ym=0;
    yp=0;
    
    
    c2=c2eq3{2}(c0,a,ym,yp,c1);
    c3=c3eq3{2}(c0,a,ym,yp,c1);
    
    hold on
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501);
        l(end+1)=plot(x,eq{3,2}(x,a(ii),c0,c1(ii),c2(ii),c3(ii)));
        l(end).DisplayName=['- c1:',num2str(c1(ii),'%.1e')];
        l(end).Color=l(ii).Color;
        l(end).LineStyle='--';
    end
    
    legend(l)
    
    % Normalised responses
    
    c2=c2eq3{1}(c0,a,ym,yp,c1);
    c3=c3eq3{1}(c0,a,ym,yp,c1);
    subplot(1,2,2)
    hold on
    clear l
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501);
        y=eq{3,1}(x,a(ii),c0,c1(ii),c2(ii),c3(ii));
        l(ii)=plot(x,y/max(abs(y)));
        l(ii).DisplayName=['+ c1:',num2str(c1(ii),'%.1e')];
    end
    
    % analysis of 3rd order (bis)
    %
    
    c2=c2eq3{2}(c0,a,ym,yp,c1);
    c3=c3eq3{2}(c0,a,ym,yp,c1);
    hold on
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501);
        y=eq{3,2}(x,a(ii),c0,c1(ii),c2(ii),c3(ii));
        l(end+1)=plot(x,y/max(abs(y)));
        l(end).DisplayName=['- c1:',num2str(c1(ii),'%.1e')];
        l(end).Color=l(ii).Color;
        l(end).LineStyle='--';
    end
    legend(l)
    
    %% (continued) Change in c1
    
    
    ym=0;
    yp=0;
    c0=0;
    c1=linspace(0,-20,7);
    a=ones(size(c1));
    
    eq3=@(x,a,c0,c1,ym,yp,eq,c2eq,c3eq) eq(x,a,c0,c1,c2eq(c0,a,ym,yp,c1),c3eq(c0,a,ym,yp,c1));
    
    figure('Name','Small response of the 3rd derivative minimisation')
    subplot(1,2,1)
    ylabel('Y Response')
    hold on
    clear l
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501);
        
        y=eq3(x,a(ii),c0,c1(ii),ym,yp,eq{3,1},c2eq3{1},c3eq3{1})...
            -eq3(x,a(ii),c0,c1(ii)+1e-5,ym,yp,eq{3,1},c2eq3{1},c3eq3{1});
        
        l(ii)=plot(x,y/max(abs(y)));
        
        l(ii).DisplayName=['+ c1:',num2str(c1(ii),'%.1e')];
    end
    
    % analysis of 3rd order negative
    
    hold on
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501);
        y=eq3(x,a(ii),c0,c1(ii),ym,yp,eq{3,2},c2eq3{2},c3eq3{2})...
            -eq3(x,a(ii),c0,c1(ii)+1e-5,ym,yp,eq{3,2},c2eq3{2},c3eq3{2});
        l(end+1)=plot(x,y/max(abs(y)));
        
        l(end).DisplayName=['- c1:',num2str(c1(ii),'%.1e')];
        l(end).Color=l(ii).Color;
        l(end).LineStyle='--';
    end
    legend(l)
    
    % Change in c1 Normal
    
    subplot(1,2,2)
    ylabel('Normal Response')
    hold on
    clear l
    kk=1;
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501)';
        
        y1=eq3(x,a(ii),c0,c1(ii),ym,yp,eq{3,1},c2eq3{1},c3eq3{1});
        y2=eq3(x,a(ii),c0,c1(ii)+1e-5,ym,yp,eq{3,1},c2eq3{1},c3eq3{1});
        [~,areaDistrib{kk}]=NormalDistance([x,y1],[x,y2]);
        l(ii)=plot(areaDistrib{kk}(:,1),areaDistrib{kk}(:,2)/max(abs(areaDistrib{kk}(:,2))));
        kk=kk+1;
        l(ii).DisplayName=['+ c1:',num2str(c1(ii),'%.1e')];
    end
    
    % analysis of 3rd order negative
    
    hold on
    for ii=1:numel(c1)
        x=linspace(-a(ii)+c0,a(ii)+c0,501)';
        y1=eq3(x,a(ii),c0,c1(ii),ym,yp,eq{3,2},c2eq3{2},c3eq3{2});
        y2=eq3(x,a(ii),c0,c1(ii)+1e-5,ym,yp,eq{3,2},c2eq3{2},c3eq3{2});
        [~,areaDistrib{kk}]=NormalDistance([x,y1],[x,y2]);
        l(end+1)=plot(areaDistrib{kk}(:,1),areaDistrib{kk}(:,2)/max(abs(areaDistrib{kk}(:,2))));
        kk=kk+1;
        
        l(end).DisplayName=['- c1:',num2str(c1(ii),'%.1e')];
        l(end).Color=l(ii).Color;
        l(end).LineStyle='--';
    end
    
    legend(l)
end


function []=HigherDerivativeOrders2()
    %% Test higher orders
    
    % Deriv1
    eq{1,1}=@(x,a,c0,c1) sqrt(a.^2-(c0.*a-x).^2)+c1;
    eq{1,2}=@(x,a,c0,c1) -sqrt(a.^2-(c0.*a-x).^2)+c1;
    
    c1eq{1,1}=@(c0,a,ym,yp) (ym-(yp-pi*a.^2/2-ym./(c0+a))./(1-1./(c0+a)))./(c0-a);
    c2eq{1,1}=@(c0,a,ym,yp) (0);
    c1eq{1,2}=@(c0,a,ym,yp) (ym-(yp+pi*a.^2/2-ym./(c0+a))./(1-1./(c0+a)))./(c0-a);
    c2eq{1,2}=@(c0,a,ym,yp) (0);
    
    % Deriv2 
    eq{2,1}=@(x,a,c0,c1,c2) a.^2./(2).*(acos(1./a.*(c0-x))...
        -(1./a.*(c0-x)).*sqrt(1-(1./a.*(c0-x)).^2))+c1.*x+c2;
    eq{2,2}=@(x,a,c0,c1,c2) -a.^2./(2).*(acos(1./a.*(c0-x))...
        -(1./a.*(c0-x)).*sqrt(1-(1./a.*(c0-x)).^2))+c1.*x+c2;
    
    c1eq{2,1}=@(c0,a,xm,xp,ym,yp,eq) (yp-ym-(eq(xp,a,c0,0,0)-eq(xm,a,c0,0,0)))./(xp-xm);
    c2eq{2,1}=@(c0,a,xm,xp,ym,yp,eq) ym-eq(xm,a,c0,0,0)-xm*(yp-ym-(eq(xp,a,c0,0,0)-eq(xm,a,c0,0,0)))./(xp-xm);
    c1eq{2,2}=c1eq{2,1};
    c2eq{2,2}=c2eq{2,1};
    
    [c3eq{1:2,1:2}]=deal(@(c0,a,ym,yp,c1) 0);
    eq2=@(a,c0,xm,xp,ym,yp,eq,c1eq,c2eq) [linspace(xm,xp,501);eq(linspace(xm,xp,501),a,c0,...
        c1eq(c0,a,xm,xp,ym,yp,eq),c2eq(c0,a,xm,xp,ym,yp,eq))];
    % Deriv 3
    eq{3,1}=@(x,a,c0,c1,c2,c3)  -a.^3./(6).*(sqrt(1-(1./a.*(c0-x)).^2).^3+sqrt(1-(1./a.*(c0-x)).^2)...
        -acos(1./a.*(c0-x)).*(1./a.*(c0-x)))+1/2*c1.*x.^2+c2.*x+c3;
    eq{3,2}=@(x,a,c0,c1,c2,c3)  -1*-a.^3./(6).*(sqrt(1-(1./a.*(c0-x)).^2).^3+sqrt(1-(1./a.*(c0-x)).^2)...
        -acos(1./a.*(c0-x)).*(1./a.*(c0-x)))+1/2*c1.*x.^2+c2.*x+c3;
    
    % Explicitely solved
%     c2eq{3,1}=@(c0,a,ym,yp,c1) (yp-ym+c1/2.*(-(c0+a).^2+(c0-a).^2)+a.^3*pi/6)./(2*a);
%     c3eq{3,1}=@(c0,a,ym,yp,c1) ym-c1/2.*(c0-a).^2-(c0-a).*(yp-ym+c1/2.*(-(c0+a).^2+(c0-a).^2)+a.^3*pi/6)./(2*a);
%     c2eq{3,2}=@(c0,a,ym,yp,c1) (yp-ym+c1/2.*(-(c0+a).^2+(c0-a).^2)-a.^3*pi/6)./(2*a);
%     c3eq{3,2}=@(c0,a,ym,yp,c1) ym-c1/2.*(c0-a).^2-(c0-a).*(yp-ym+c1/2.*(-(c0+a).^2+(c0-a).^2)-a.^3*pi/6)./(2*a);
    

    c2eq{3,1}=@(c0,a,xm,xp,ym,yp,c1,eq) (yp-ym-(eq(xp,a,c0,c1,0,0)-eq(xm,a,c0,c1,0,0)))./(xp-xm);
    c3eq{3,1}=@(c0,a,xm,xp,ym,yp,c1,eq) ym-eq(xm,a,c0,c1,0,0)-xm*(yp-ym-(eq(xp,a,c0,c1,0,0)-eq(xm,a,c0,c1,0,0)))./(xp-xm);
    c2eq{3,2}=c2eq{3,1};
    c3eq{3,2}=c3eq{3,1};

    %eq3=@(x,a,c0,c1,ym,yp,eq,c2eq,c3eq) eq(x,a,c0,c1,c2eq(c0,a,ym,yp,c1),c3eq(c0,a,ym,yp,c1));
    eq3=@(a,c0,c1,xm,xp,ym,yp,eq,c2eq,c3eq) [linspace(xm,xp,501);eq(linspace(xm,xp,501),a,c0,...
        c1,c2eq(c0,a,xm,xp,ym,yp,c1,eq),c3eq(c0,a,xm,xp,ym,yp,c1,eq))];
    
    plotFuncs{1}= @(p,pd,pn) plot(p(:,1),p(:,2));
    plotFuncs{2}= @(p,pd,pn) plot(p(:,1),p(:,2)/max(abs(p(:,2))));
    plotFuncs{3}= @(p,pd,pn) plot(p(:,1),(p(:,2)-pd(:,2))/max(abs(pd(:,2)-p(:,2))));
    plotFuncs{4}= @(p,pd,pn) plot(pn(:,1),pn(:,2)/max(abs(pn(:,2))));
    %% Analysis of 3rd Order using oneshot func
    ii=3;
    cellarg3=[eq(ii,:)',c2eq(ii,:)',c3eq(ii,:)'];
    nameCell={'+','-'};
    nameCol=2;
    
    
    c0=linspace(0,1-1e-5,7);
    a=2;
    c1=0;
    xm=-1;
    xp=1;
    ym=0;
    yp=0;
    
    numarg=BuildFullFactorialMatrix(a,c0,c1,xm,xp,ym,yp);
    %ExploreFunction(eq3,numarg,cellarg3,plotFuncs,nameCol,nameCell)
    
   %% Analysis of 2nd Order using oneshot func
    ii=2;
    cellarg2=[eq(ii,:)',c1eq(ii,:)',c2eq(ii,:)'];
    nameCell={'+','-'};
    nameCol=2;
    
    
    xm=-1;
    xp=1;
    ym=0;
    yp=0;
    
    c0=linspace(-1,1-01e-5,7);
    a=2;
    
    
    numarg=BuildFullFactorialMatrix(a,c0,xm,xp,ym,yp);
    %numarg(:,2)=numarg(:,1)-xp-1e-5;
    ExploreFunction(eq2,numarg,cellarg2,plotFuncs,nameCol,nameCell)
end