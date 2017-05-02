function [matCoeff,pp]=PiecePolyBasis()
    
    d=1;
    Dv=1.68;
    steps=15;
    plotSteps=20;
    precision=8;
    
    [matCoeff]=GenerateMatCoeff(d,Dv,steps,precision);
    
        PlotPiecePoly(d,Dv,plotSteps,matCoeff);
        figure,
        semilogy(1:steps,abs(matCoeff))
    
    [pp]=CalculatePolyObj(matCoeff,d);
end

function [pp]=CalculatePolyObj(matCoeff,d)
    
    calcMask=@(e) [1 0 0;-2*e 1 0; e^2 -e 1];
    
    nStep=size(matCoeff,2);
    realCoeff=zeros(size(matCoeff));
    
    %     for ii=0:nStep-1
    %         e=(2*ii-1)*d;
    %         realCoeff(:,ii+1)=calcMask(e)*matCoeff(:,ii+1);
    %     end
    negCoeff=zeros(size(matCoeff(:,2:end)));
    for ii=nStep:-1:2
        negCoeff(:,nStep+1-ii)=calcMask(2*d)*([1;-1;1].*matCoeff(:,ii));
    end
    
    ppBounds=(-2*d*(nStep-1)-d):2*d:(2*d*(nStep-1)+d);
    allCoeff=[negCoeff,matCoeff];
    allCoeff=double(allCoeff);
    pp=mkpp(ppBounds,allCoeff');
    x=linspace(-99,99,10000);
    figure, plot(x,ppval(x,pp))
end

function []=PlotPiecePoly(d,Dv,plotSteps,matCoeff)
    
    x=linspace(0,2*d,plotSteps+1);
    x(end)=[];
    
    R=[x'.^2,x',ones(size(x'))];
    Rderiv=[2*x',ones(size(x')),zeros(size(x'))];
    Rintg=[x'.^3/3,x'.^2/2,x'];
    
    Y=R*matCoeff;
    Yderiv=Rderiv*matCoeff;
    Yintg=Rintg*matCoeff;
    Yintg(:,2:end)=Yintg(:,2:end)+Dv;
    
    [Yvec]=VectoriseMatrix(Y);
    [YderivV]=VectoriseMatrix(Yderiv);
    [YintegV]=VectoriseMatrix(Yintg);
    
    Xvec=zeros([1,numel(Y)]);
    colLength=size(Y,1);
    for ii=1:size(Y,2)
        iS=1+colLength*(ii-1);
        iE=colLength*ii;
        Xvec(iS:iE)=x+((-1+2*(ii-1))*d);
    end
    Xvec=[-Xvec(end:-1:colLength),Xvec];
    YderivV=[-YderivV(end:-1:colLength),YderivV];
    Yvec=[Yvec(end:-1:colLength),Yvec];
    YintegV=[-YintegV(end:-1:colLength)+Dv,YintegV];
    
    
    figure
    plot(Xvec,Yvec,Xvec,YderivV,Xvec,YintegV)
    grid on
    
end

function [Yvec]=VectoriseMatrix(Y)
    
    Yvec=zeros([1,numel(Y)]);
    colLength=size(Y,1);
    for ii=1:size(Y,2)
        iS=1+colLength*(ii-1);
        iE=colLength*ii;
        Yvec(iS:iE)=Y(:,ii)';
    end
    
end

function [matCoeff]=GenerateMatCoeff(d,Dv,steps,precision)
    
    matCoeff=vpa(zeros([3,steps]),precision);
    %matCoeff=(zeros([3,steps]));
    [matCoeff(1,1),matCoeff(2,1),matCoeff(3,1)]=GetOriginalCoeff(d,Dv,precision);
    
    for ii=2:steps
        [matCoeff(1,ii),matCoeff(2,ii),matCoeff(3,ii)]=...
            NextCoeff(matCoeff(1,ii-1),matCoeff(2,ii-1),matCoeff(3,ii-1),d);
    end
    
end

function [an,bn,cn]=NextCoeff(ai,bi,ci,d)
    
    
    
    cn=4*ai*d^2+2*bi*d+ci;
    bn=4*ai*d+bi;
    an=3*(-(2*bn*d^2+2*cn*d))/(8*d^3);
    %an=3/(8*d^3)*((8/3*d^3*ai+2*bi*d^2+2*ci*d)-(2*bi*d^2+2*ci*d));
    %    (8/3*d^3*ai+2*bi*d^2+2*ci*d)
    %    (8/3*d^3*an+2*bn*d^2+2*cn*d)
    
end

function [coeffMask,endCond,condVec]=MatrixOpCoeff(d,precision)
    
    coeffMask=[-6, -9/(4*d), -3/(4*d^2);
        4*d, 1,  0
        4*d^2, 2*d, 1];
    endCond=[4*d^2 2*d  1];
    [V,l]=eig(vpa(coeffMask,precision));
    %[V,l]=eig(coeffMask);
    [~,x]=max(max(abs(l)));
    eigLim=zeros(3);
    eigLim(x,x)=1;
    
    condVec=endCond*V*eigLim/V;
    %     for ii=1:10,
    %         res=endCond*(coeffMask^ii);
    %         strRes(ii,1:3)=res/res(end);
    %     end
    %     figure
    %     plot(1:10,strRes')
end

function [a0,b0,c0]=GetOriginalCoeff(d,Dv,precision)
    
    
    
    [~,~,l]=MatrixOpCoeff(d,precision);
    %a0=(Dv-2*c0*d)/(8/3*d^3-4*d^3);
    
    a0=-Dv*l(3)/(2*d)/(l(1)-2*d*l(2)+2/3*d^2*l(3));
    
    c0=+2/3*d^2*a0+Dv/(2*d);
    b0=-2*a0*d;
    
end

function [lineSum]=PlotBasisSums(desRange,d,pp,basisInd,coeffs)
    
    if nargin<5; coeffs=ones(size(basisInd));end
    
    x=linspace(desRange(1),desRange(2),10001);
    
    kk=1;
    for ii=basisInd;
        YMat(kk,1:numel(x))=ppval(pp,x+(ii)*2*d)*coeffs(kk);
        kk=kk+1;
    end
    figure
    plot(x,YMat)
    hold on, 
    lineSum=plot(x,sum(YMat),'r--');
    
end