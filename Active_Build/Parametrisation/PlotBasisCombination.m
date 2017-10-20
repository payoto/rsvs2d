

%% Coeff Ratios
nBases=9;
matCoeff=double(matCoeff);
desRange=[-50 50];basisInd=-nBases:nBases;

coeffs=matCoeff(1,1:nBases+1)/matCoeff(1,1);
coeffs(2:end)=-coeffs(2:end);
coeffs=coeffs([end:-1:1,2:end]);
PlotBasisSums(desRange,d,pp,basisInd,coeffs)

%% Parallel Peak Scaling

nBases=1;
matCoeff=double(matCoeff);
desRange=[-50 50];basisInd=-nBases:nBases;

xpeak=@(a,b) -b./(2*a);
peakHeight=@(a,b,c,xp) a.*xp.^2+b.*xp+c;
a=matCoeff(1,1:nBases+1);
b=matCoeff(2,1:nBases+1);
c=matCoeff(3,1:nBases+1);

coeffs=peakHeight(a,b,c,xpeak(a,b));
coeffs=coeffs/coeffs(1);
coeffs(2:end)=-coeffs(2:end);
coeffs=coeffs([end:-1:1,2:end]);
PlotBasisSums(desRange,d,pp,basisInd,coeffs)

%% Sequential Peak Scaling

nBases=10;
matCoeff=double(matCoeff);

desRange=[-30 30];
xpeak=@(a,b) -b./(2*a);
peakHeight=@(a,b,c,xp) a.*xp.^2+b.*xp+c;
basisInd=0;
coeffs=1;
basePeak=peakHeight(matCoeff(1,1),matCoeff(2,1),matCoeff(3,1),...
xpeak(matCoeff(1,1),matCoeff(2,1)));
legEntry={};

lBasis(1)=PlotBasisSums(desRange,d,pp,-1:1,[0 1 0]);
legEntry{1}=['Smoothing Level: ',int2str(0)];
for ii=1:nBases
    
    actPoly=ii-basisInd+1;
    a=sum(matCoeff(1,actPoly).*coeffs);
    b=sum(matCoeff(2,actPoly).*coeffs);
    c=sum(matCoeff(3,actPoly).*coeffs);
    xp=xpeak(a,b);
    peakHeightVal=peakHeight(a,b,c,xpeak(a,b));
    offSetRatio=(d+abs(d-xp))/d;
    
    newCoeff=-peakHeightVal/basePeak*offSetRatio;
    
    coeffs=[newCoeff,coeffs,newCoeff];
    
    basisInd=-ii:ii;
    lBasis(ii+1)=PlotBasisSums(desRange,d,pp,basisInd,coeffs);
    
    legEntry{ii+1}=['Smoothing Level: ',int2str(ii)];
end
h=figure('Name','Sequential Peak basis smoothing');
ax=subplot(1,2,1);
lNew=copyobj(lBasis,ax);
cOrd=get(gca,'ColorOrder');
for ii=1:length(lNew)
    lNew(ii).LineStyle='-';
    lNew(ii).Color=cOrd(mod(ii-1,7)+1,:);
end
title('Smoothed basis function')
legend(lNew,legEntry)
overShoot=[];
for ii=1:length(lNew)
    overShoot(ii)=abs(min(lNew(ii).YData))/abs(max(lNew(ii).YData));
end
subplot(2,2,4)
semilogy(0:nBases,overShoot)
xlabel('Smoothing Level')
ylabel('Overshoot/peak')
title('Overshoot Comparison')
subplot(2,2,2)
hold on
lGrad=[];
for ii=1:length(lNew)
    
    grad=[lNew(ii).YData(1)-lNew(ii).YData(2),...
        lNew(ii).YData(1:end-2)-lNew(ii).YData(3:end),...
        lNew(ii).YData(end-1)-lNew(ii).YData(end)]./...
        [lNew(ii).XData(1)-lNew(ii).XData(2),...
        lNew(ii).XData(1:end-2)-lNew(ii).XData(3:end),...
        lNew(ii).XData(end-1)-lNew(ii).XData(end)];
    lGrad(ii)=plot(lNew(ii).XData,grad);
end
title('Derivatives')
legend(lGrad,legEntry)
coeffstruct.peaksmooth.coeff=coeffs;
coeffstruct.peaksmooth.maxsmooth=10;
%% Sequential Peak Scaling 2

nBases=10;
matCoeff=double(matCoeff);

desRange=[-30 30];
xpeak=@(a,b) -b./(2*a);
peakHeight=@(a,b,c,xp) a.*xp.^2+b.*xp+c;
basisInd=0;
coeffs=1;
basePeak=peakHeight(matCoeff(1,1),matCoeff(2,1),matCoeff(3,1),...
xpeak(matCoeff(1,1),matCoeff(2,1)));
legEntry={};

lBasis(1)=PlotBasisSums(desRange,d,pp,-1:1,[0 1 0]);
legEntry{1}=['Smoothing Level: ',int2str(0)];
for ii=1:nBases
    
    actPoly=ii-basisInd+1;
    a=sum(matCoeff(1,actPoly).*coeffs);
    b=sum(matCoeff(2,actPoly).*coeffs);
    c=sum(matCoeff(3,actPoly).*coeffs);
    xp=xpeak(a,b);
    peakHeightVal=peakHeight(a,b,c,xpeak(a,b));
    bCoeff=matCoeff(:,1)+matCoeff(:,1+ii*2);
    basePeak=peakHeight(bCoeff(1,1),bCoeff(2,1),bCoeff(3,1),...
        xpeak(bCoeff(1,1),bCoeff(2,1)));
    bxp=xpeak(bCoeff(1),bCoeff(2));
    offSetRatio=(bxp+abs(bxp-xp))/(bxp+abs(bxp-d));
    
    newCoeff=-peakHeightVal/basePeak*offSetRatio;
    
    coeffs=[newCoeff,coeffs,newCoeff];
    
    basisInd=-ii:ii;
    lBasis(ii+1)=PlotBasisSums(desRange,d,pp,basisInd,coeffs);
    
    legEntry{ii+1}=['Smoothing Level: ',int2str(ii)];
end
h=figure('Name','Sequential Peak basis smoothing');
ax=subplot(1,2,1);
lNew=copyobj(lBasis,ax);
cOrd=get(gca,'ColorOrder');
for ii=1:length(lNew)
    lNew(ii).LineStyle='-';
    lNew(ii).Color=cOrd(mod(ii-1,7)+1,:);
end
title('Smoothed basis function')
legend(lNew,legEntry)
overShoot=[];
for ii=1:length(lNew)
    overShoot(ii)=abs(min(lNew(ii).YData))/abs(max(lNew(ii).YData));
end
subplot(2,2,4)
semilogy(0:nBases,overShoot)
xlabel('Smoothing Level')
ylabel('Overshoot/peak')
title('Overshoot Comparison')
subplot(2,2,2)
hold on
lGrad=[];
for ii=1:length(lNew)
    
    grad=[lNew(ii).YData(1)-lNew(ii).YData(2),...
        lNew(ii).YData(1:end-2)-lNew(ii).YData(3:end),...
        lNew(ii).YData(end-1)-lNew(ii).YData(end)]./...
        [lNew(ii).XData(1)-lNew(ii).XData(2),...
        lNew(ii).XData(1:end-2)-lNew(ii).XData(3:end),...
        lNew(ii).XData(end-1)-lNew(ii).XData(end)];
    lGrad(ii)=plot(lNew(ii).XData,grad);
end
title('Derivatives')
legend(lGrad,legEntry)
coeffstruct.peaksmooth.coeff=coeffs;
coeffstruct.peaksmooth.maxsmooth=10;

%% Sequential Peak Scaling 3

nBases=10;
matCoeff=double(matCoeff);

desRange=[-30 30];
xpeak=@(a,b) -b./(2*a);
peakHeight=@(a,b,c,xp) a.*xp.^2+b.*xp+c;
basisInd=0;
coeffs=1;
basePeak=peakHeight(matCoeff(1,1),matCoeff(2,1),matCoeff(3,1),...
xpeak(matCoeff(1,1),matCoeff(2,1)));
legEntry={};

lBasis(1)=PlotBasisSums(desRange,d,pp,-1:1,[0 1 0]);
legEntry{1}=['Smoothing Level: ',int2str(0)];
for ii=1:nBases
    
    actPoly=ii-basisInd+1;
    actPoly=actPoly+1;
    a=sum(matCoeff(1,actPoly).*coeffs);
    b=sum(matCoeff(2,actPoly).*coeffs);
    c=sum(matCoeff(3,actPoly).*coeffs);
    xp=xpeak(a,b);
    peakHeightVal=peakHeight(a,b,c,xpeak(a,b));
    
    currCoeff=matCoeff(:,2)+matCoeff(:,2+ii*2);
    
    bCoeff=matCoeff(:,1)+matCoeff(:,1+ii*2);
    basePeak=peakHeight(bCoeff(1,1),bCoeff(2,1),bCoeff(3,1),...
        xpeak(bCoeff(1,1),bCoeff(2,1)));
    bxp=xpeak(bCoeff(1),bCoeff(2));
    offSetRatio=(bxp+abs(bxp-xp))/(bxp+abs(bxp-d));
    
    newCoeff=-peakHeightVal/basePeak*offSetRatio;
    
    newCoeff=-mean([a,b,c]'./currCoeff);
    
    coeffs=[newCoeff,coeffs,newCoeff];
    
    basisInd=-ii:ii;
    lBasis(ii+1)=PlotBasisSums(desRange,d,pp,basisInd,coeffs);
    
    legEntry{ii+1}=['Smoothing Level: ',int2str(ii)];
end
h=figure('Name','Sequential Peak basis smoothing');
ax=subplot(1,2,1);
lNew=copyobj(lBasis,ax);
cOrd=get(gca,'ColorOrder');
for ii=1:length(lNew)
    lNew(ii).LineStyle='-';
    lNew(ii).Color=cOrd(mod(ii-1,7)+1,:);
end
title('Smoothed basis function')
legend(lNew,legEntry)
overShoot=[];
for ii=1:length(lNew)
    overShoot(ii)=abs(min(lNew(ii).YData))/abs(max(lNew(ii).YData));
end
subplot(2,2,4)
semilogy(0:nBases,overShoot)
xlabel('Smoothing Level')
ylabel('Overshoot/peak')
title('Overshoot Comparison')
subplot(2,2,2)
hold on
lGrad=[];
for ii=1:length(lNew)
    
    grad=[lNew(ii).YData(1)-lNew(ii).YData(2),...
        lNew(ii).YData(1:end-2)-lNew(ii).YData(3:end),...
        lNew(ii).YData(end-1)-lNew(ii).YData(end)]./...
        [lNew(ii).XData(1)-lNew(ii).XData(2),...
        lNew(ii).XData(1:end-2)-lNew(ii).XData(3:end),...
        lNew(ii).XData(end-1)-lNew(ii).XData(end)];
    lGrad(ii)=plot(lNew(ii).XData,grad);
end
title('Derivatives')
legend(lGrad,legEntry)
coeffstruct.peaksmooth.coeff=coeffs;
coeffstruct.peaksmooth.maxsmooth=10;

%% Poly coeff smoothing
matCoeff=double(matCoeff);
nBases=30;
%nSmooth=[0 1 2 3 4 5 6 7 8 9 10];
nSmooth=[0 1 2 3 4 5 6 7 8];
desRange=[-50 50];
basisInd=-nBases:nBases;
lBasis=[];
legEntry={};
for kk=1:length(nSmooth)
    if nSmooth(kk)>0
        coeffs=matCoeff(1,1:nBases+1)/matCoeff(1,1);
        coeffs(2:end)=-coeffs(2:end);
        coeffs=coeffs([end:-1:1,2:end]);
        
        for lll=2:nSmooth(kk);
            coeffMat2=zeros((length(coeffs)*[1 1]+[1 0]));
            coeffMat2(1,:)=coeffs;
            for ii=1:length(coeffs)
                
                iSD=max([ii-nBases,1]);
                iED=min([ii+nBases,length(coeffs)]);
                iSS=max([nBases+2-ii,1]);
                iES=iSS+(iED-iSD);
                
                coeffMat2(ii+1,iSD:iED)=coeffs(ii)*coeffs(iSS:iES);
                
            end
            
            coeffs=sum(coeffMat2(2:end,:));
            coeffs=coeffs/coeffs(nBases+1);
        end
    else
        coeffs=zeros([1,2*nBases+1]);
        coeffs(nBases+1)=1;
    end
    coeffs(abs(coeffs)<1e-3)=0;
    coeffstruct.polysmooth(kk).coeff=coeffs(coeffs~=0);
    coeffstruct.polysmooth(kk).maxsmooth=nSmooth(kk);
    lBasis(kk)=PlotBasisSums(desRange,d,pp,basisInd,coeffs);
    legEntry{kk}=['Smoothing Level: ',int2str(nSmooth(kk))];
end

h=figure('Name','Polynomial Coefficient based smoothing');
ax=subplot(1,2,1);
lNew=copyobj(lBasis,ax);
cOrd=get(gca,'ColorOrder');
for ii=1:length(lNew)
    lNew(ii).LineStyle='-';
    lNew(ii).Color=cOrd(mod(ii-1,7)+1,:);
end
title('Smoothed basis function')
legend(lNew,legEntry)
overShoot=[];
for ii=1:length(lNew)
    overShoot(ii)=abs(min(lNew(ii).YData))/abs(max(lNew(ii).YData));
end
subplot(2,2,4)
semilogy(nSmooth,overShoot)
xlabel('Smoothing Level')
ylabel('Overshoot/peak')
title('Overshoot Comparison')
subplot(2,2,2)
hold on
lGrad=[];
for ii=1:length(lNew)
    
    grad=[lNew(ii).YData(1)-lNew(ii).YData(2),...
        lNew(ii).YData(1:end-2)-lNew(ii).YData(3:end),...
        lNew(ii).YData(end-1)-lNew(ii).YData(end)]./...
        [lNew(ii).XData(1)-lNew(ii).XData(2),...
        lNew(ii).XData(1:end-2)-lNew(ii).XData(3:end),...
        lNew(ii).XData(end-1)-lNew(ii).XData(end)];
    lGrad(ii)=plot(lNew(ii).XData,grad);
end
title('Derivatives')
legend(lGrad,legEntry)
%% Poly peak smoothing
matCoeff=double(matCoeff);
nBases=30;
%nSmooth=[0 1 2 3 4 5 6 7 8 9 10];
nSmooth=[0 1 2 3 4 5 6 7 8];
desRange=[-50 50];
basisInd=-nBases:nBases;
lBasis=[];
legEntry={};
for kk=1:length(nSmooth)
    if nSmooth(kk)>0
        xpeak=@(a,b) -b./(2*a);
        peakHeight=@(a,b,c,xp) a.*xp.^2+b.*xp+c;
        a=matCoeff(1,1:nBases+1);
        b=matCoeff(2,1:nBases+1);
        c=matCoeff(3,1:nBases+1);
        
        coeffs=peakHeight(a,b,c,xpeak(a,b));
        coeffs=coeffs/coeffs(1);
        coeffs(2:end)=-coeffs(2:end);
        coeffs=coeffs([end:-1:1,2:end]);
        
        for lll=2:nSmooth(kk);
            coeffMat2=zeros((length(coeffs)*[1 1]+[1 0]));
            coeffMat2(1,:)=coeffs;
            for ii=1:length(coeffs)
                
                iSD=max([ii-nBases,1]);
                iED=min([ii+nBases,length(coeffs)]);
                iSS=max([nBases+2-ii,1]);
                iES=iSS+(iED-iSD);
                
                coeffMat2(ii+1,iSD:iED)=coeffs(ii)*coeffs(iSS:iES);
                
            end
            
            coeffs=sum(coeffMat2(2:end,:));
            coeffs=coeffs/coeffs(nBases+1);
        end
    else
        coeffs=zeros([1,2*nBases+1]);
        coeffs(nBases+1)=1;
    end
    coeffs(abs(coeffs)<1e-3)=0;
    
    coeffstruct.polypeaksmooth(kk).coeff=coeffs(coeffs~=0);
    coeffstruct.polypeaksmooth(kk).maxsmooth=nSmooth(kk);
    lBasis(kk)=PlotBasisSums(desRange,d,pp,basisInd,coeffs);
    legEntry{kk}=['Smoothing Level: ',int2str(nSmooth(kk))];
end

h=figure('Name','Parrallel Peak basis smoothing');
ax=subplot(1,2,1);
lNew=copyobj(lBasis,ax);
cOrd=get(gca,'ColorOrder');
for ii=1:length(lNew)
    lNew(ii).LineStyle='-';
    lNew(ii).Color=cOrd(mod(ii-1,7)+1,:);
end
title('Smoothed basis function')
legend(lNew,legEntry)
overShoot=[];
for ii=1:length(lNew)
    overShoot(ii)=abs(min(lNew(ii).YData))/abs(max(lNew(ii).YData));
end
subplot(2,2,4)
semilogy(nSmooth,overShoot)
xlabel('Smoothing Level')
ylabel('Overshoot/peak')
title('Overshoot Comparison')
subplot(2,2,2)
hold on
lGrad=[];
for ii=1:length(lNew)
    
    grad=[lNew(ii).YData(1)-lNew(ii).YData(2),...
        lNew(ii).YData(1:end-2)-lNew(ii).YData(3:end),...
        lNew(ii).YData(end-1)-lNew(ii).YData(end)]./...
        [lNew(ii).XData(1)-lNew(ii).XData(2),...
        lNew(ii).XData(1:end-2)-lNew(ii).XData(3:end),...
        lNew(ii).XData(end-1)-lNew(ii).XData(end)];
    lGrad(ii)=plot(lNew(ii).XData,grad);
end
title('Derivatives')
legend(lGrad,legEntry)
%% Comparison to normal distribution
if false %
    figure
    subplot(1,2,1)
    hold on
    for ii=1:length(lNew)
        
        normApprox=max(lNew(ii).YData)*pdf('Normal',lNew(ii).XData,0,nSmooth(ii)*2+1)...
            /(pdf('Normal',0,0,nSmooth(ii)*2+1));
        
        err(ii,1:numel(lNew(ii).YData))=(lNew(ii).YData...
            -normApprox);
        lErr(ii)=plot(lNew(ii).XData,err(ii,1:numel(lNew(ii).YData)));
        lErr(ii).Color=lNew(ii).Color;
        lErr(ii).DisplayName=lNew(ii).DisplayName;
        
        rmsErr(ii)=sqrt(sum((lNew(ii).YData-normApprox).^2))/numel(lNew(ii).YData);
    end
    legend(lErr,'Location','SouthEast')
    subplot(1,2,2)
    plot(nSmooth,rmsErr)
end

%% Study arc length

ParaBolaArcLength=@(a,b,x) (sqrt((2*a*x+b).^2+1).*(2*a*x+b)+asinh(2*a*x+b))/(4*a);

desRange=[-10 10];

x=linspace(desRange(1),desRange(2),5001);
