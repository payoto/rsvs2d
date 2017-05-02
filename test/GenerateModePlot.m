%%

integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*(tDistrib(1:end-1)+tDistrib(2:end))/2]);
keepPos=find(l2.YData<-0.017);

xDat=l(1).XData;
yDat=vertcat(l(:).YData);
%%
[~,xPos]=max(yDat(1,:));
xDat=xDat/(xDat(xPos+2)-xDat(xPos));
%%
[~,xPos]=max(yDat(1,:));
xDat=xDat-xDat(xPos);

yDat(1,:)=yDat(1,:)/max(max(yDat(1,:)));
yDat(2,:)=yDat(2,:)/max(max(yDat(2,:)));
yDat(3,:)=[];
yDeriv=([yDat(:,2:end),[0;0]]-[[0;0],yDat(:,1:end-1)])./repmat(([xDat(2:end),0]-[0,xDat(1:end-1)]),[2,1]);
yInteg=integr(xDat,yDat(1,:));
yInteg(2,:)=integr(xDat,yDat(2,:));
yDat(3,:)=yDat(2,:)-yDat(1,:);
yInteg(3,:)=yInteg(2,:)-yInteg(1,:);
yDeriv(3,:)=yDeriv(2,:)-yDeriv(1,:);
for ii=1:3;
    figure,
    hold on,
    plot(xDat(1,keepPos),yDat(ii,keepPos))
    plot(xDat(1,keepPos),yDeriv(ii,keepPos))
    plot(xDat(1,keepPos),yInteg(ii,keepPos))
end


%%
figNames={'Snake2','Analytical','Difference Anal-Snak'}

for ii=1:3;
    h=findobj(ii);
    figure(ii)
    h(1).PaperPositionMode='auto';
    h(1).Name=figNames{ii};
    h(2).XLim=[-20 20];
    h(1).Position=[200 200 400 300];
    plot([-1 -1],[-1 2],'Color',[1 1 1]*0.3,'LineStyle','--')
    plot([1 1],[-1 2],'Color',[1 1 1]*0.3,'LineStyle','--')
    l=findobj(ii,'type','line');
    plot([-20,-0,0,20],[0,0,l(3).YData(end),l(3).YData(end)],'Color',l(3).Color,'LineStyle','--')
end


for ii=1:3;
    h=findobj(ii);
    print(h(1),['.\fig\',h(1).Name,'.eps'],'-depsc')
end
for ii=1:3;
    h=findobj(ii);
    hgsave(h(1),['.\fig\',h(1).Name,'.fig'])
end

%%

nFig=6

for ii=1:nFig;
    h=findobj(ii,'type','axes');
    if numel(h)>0
        for jj=1:numel(h)
            h(jj).TickLabelInterpreter='latex';
            h(jj).XLabel.Interpreter='latex';
            h(jj).YLabel.Interpreter='latex';
        end
    end
    h=findobj(ii,'type','legend');
    if numel(h)>0
        for jj=1:numel(h)
            h(jj).Interpreter='latex';
        end
    end
end
%%
for ii=1:nFig;
    h=findobj(ii);
    print(h(1),['.\fig\',h(1).Name,'.eps'],'-depsc')
    hgsave(h(1),['.\fig\',h(1).Name,'.fig'])
end
