

%% Coeff Ratios
nBases=9;
matCoeff=double(matCoeff);
desRange=[-50 50];basisInd=-nBases:nBases;

coeffs=matCoeff(1,1:nBases+1)/matCoeff(1,1);
coeffs(2:end)=-coeffs(2:end);
coeffs=coeffs([end:-1:1,2:end]);
PlotBasisSums(desRange,d,pp,basisInd,coeffs)

%%
matCoeff=double(matCoeff);
nBases=8;
%nSmooth=[0 1 2 3 4 5 6 7 8 9 10];
nSmooth=[0 1 2 3 4 5 6 7 8];
desRange=[-30 30];
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
    coeffs
    lBasis(kk)=PlotBasisSums(desRange,d,pp,basisInd,coeffs);
    legEntry{kk}=['Smoothing Level: ',int2str(nSmooth(kk))];
end

h=figure('Name','Smoothed Basis');
ax=axes;
lNew=copyobj(lBasis,ax);
cOrd=get(gca,'ColorOrder');
for ii=1:length(lNew)
    lNew(ii).LineStyle='-';
    lNew(ii).Color=cOrd(mod(ii-1,7)+1,:);
end
legend(lNew,legEntry)
overShoot=[];
for ii=1:length(lNew)
    overShoot(ii)=abs(min(lNew(ii).YData))/abs(max(lNew(ii).YData));
end
figure('Name','Overshoot Comparison')
semilogy(nSmooth,overShoot)
xlabel('Smoothing Level')
ylabel('Overshoot/peak')




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


