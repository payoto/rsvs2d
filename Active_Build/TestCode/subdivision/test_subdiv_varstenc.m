function [figH2, figH, l]=test_subdiv_varstenc(points,strct, maxRef, figH2, isplot)

if ~exist('isplot','var')
    isplot = true;
end
    
d2pds2=@(p,pp,pm) (2*(-p*(norm(p-pp)+norm(p-pm))+norm(p-pm)*pp+norm(p-pp)*pm)/((norm(p-pp)*norm(p-pm)^2)+(norm(p-pp)^2*norm(p-pm))));
 normd2pds2=@(p,pp,pm) norm(d2pds2(p,pp,pm));
oldPoints=[0 0;1 0; 1 1;0 1];
if isplot
    figH = figure;
    subplot(1,2,1)
    hold on
else
    figH = 0;
end
newPoints=points;
if isplot
    plot(newPoints([1:end,1],1),newPoints([1:end,1],2),'o-')
    subplot(1,2,2)
    hold on,
end
kk=1;
vec=1:maxRef;
for ii=vec,
    [newPoints]=SubSurfVarStencil_NoCorn_STA(points,ii,strct);
    if isplot
        subplot(1,2,1)
        plot(newPoints(:,1),newPoints(:,2))
    end
    for jj=2:length(newPoints(:,1))-1;
        normCurv(jj-1)=normd2pds2(newPoints(jj,:),newPoints(jj-1,:),newPoints(jj+1,:));
    end
    if isplot
        subplot(1,2,2)
    
        h(kk)=plot((0:(length(newPoints(:,1))-3))/(length(newPoints(:,1))-3),normCurv);
    end
    cellLeg{kk}=num2str(ii);
    
    meanCurv(kk)=mean(normCurv);
    maxCurv(kk)=max(abs(normCurv));
    stdCurv(kk)=std(normCurv);
    kk=kk+1;
  
end
if isplot
    subplot(1,2,1)
    legend(h,cellLeg);
end
if ~exist('figH2','var')
    figH2 = figure;
end
figure(figH2)
subplot(1,3,1)
hold on
l(1) = plot(vec,meanCurv);
subplot(1,3,2)
hold on
l(2) = plot(vec,stdCurv);
subplot(1,3,3)
hold on
l(3) = plot(vec,maxCurv);
end