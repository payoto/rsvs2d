figure
for ll=1:10
clear h
clear cellLeg
clear meanCurv,
clear stdCurv
clear normCurv
d2pds2=@(p,pp,pm) (2*(-p*(norm(p-pp)+norm(p-pm))+norm(p-pm)*pp+norm(p-pp)*pm)/((norm(p-pp)*norm(p-pm)^2)+(norm(p-pp)^2*norm(p-pm))))
 normd2pds2=@(p,pp,pm) norm(d2pds2(p,pp,pm))
oldPoints=[0 0; 1 0; 1 1;0 1];
oldPoints=[oldPoints;[0 0;0 -1;0 -2; 0 -3]];
oldPoints=[oldPoints;[-2 -3;-2 0]];

figure
subplot(1,2,1)
hold on

newPoints=oldPoints;

subplot(1,2,2)
hold on,
kk=1;
global testCoeffTest;
vec=0.0:0.0005:0.1;
for ii=vec,
    testCoeffTest=ii;
    newPoints=SubDivision(oldPoints,ll,'homemade');
    subplot(1,2,1)
    %plot(newPoints(:,1),newPoints(:,2))
    for jj=2:length(newPoints(:,1))-1;
        normCurv(jj-1)=normd2pds2(newPoints(jj,:),newPoints(jj-1,:),newPoints(jj+1,:));
    end,
    subplot(1,2,2)
    %h(kk)=plot(1:(length(newPoints(:,1))-2),normCurv);
    cellLeg{kk}=num2str(ii);
    
    meanCurv(kk)=mean(normCurv);
    stdCurv(kk)=std(normCurv);
    kk=kk+1;
  
end
subplot(1,2,1)
%legend(h,cellLeg);
figure(1)
subplot(1,2,1)
hold on
plot(vec,meanCurv)
subplot(1,2,2)
hold on
plot(vec,stdCurv)
pause(1)
end