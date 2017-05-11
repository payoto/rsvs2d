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

for ii=figList
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


for ii=figList
    h=findobj(ii);
    print(h(1),['.\fig\',h(1).Name,'.eps'],'-depsc')
end
for ii=figList
    h=findobj(ii);
    hgsave(h(1),['.\fig\',h(1).Name,'.fig'])
end

%%


for ii=figList;
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
    h=findobj(ii,'type','text');
    if numel(h)>0
        for jj=1:numel(h)
            h(jj).Interpreter='latex';
        end
    end
end
%%
for ii=figList;
    h=findobj(ii);
    figure(ii)
    h(1).PaperPositionMode='auto';
%     h(1).Position=[100 100 450 325];
end

%%
for ii=figList;
    h=findobj(ii,'type','axes');
    if numel(h)>0
        for jj=1:numel(h)
            h(jj).Position=[0.1300 0.1100 0.7750 0.8150];
            h(jj).OuterPosition=[0 0 1 1];
        end
    end
end
%%
for ii=figList;
    h=findobj(ii);
    h(1).Renderer='painters';
    h(1).Color='none';
    print(h(1),'-r300','-dpng',['.\fig\',h(1).Name,'.png'])
    print(h(1),'-depsc',['.\fig\',h(1).Name,'.eps'])
    h(1).Color=[0.93 0.93 0.93];
    hgsave(h(1),['.\fig\',h(1).Name,'.fig'])
end

%%

for ii=figList;
    h=findobj(ii,'type','axes');
    if numel(h)>0
        for jj=1:numel(h)
            h(jj).Position=hRef(3).Position;
            h(jj).OuterPosition=hRef(3).OuterPosition;
        end
    end
    h=findobj(ii);
    h(1).Position=hRef(1).Position;
end

%%
for ii=figList
    
    l=findobj(ii,'type','line');
    yDat=vertcat(l(:).YData);
    maxDat=max(-yDat,[],2);
    maxDat=maxDat(maxDat>0);
    sumdat(ii).dat=maxDat;
    sumdat(ii).mean=mean(maxDat);
    sumdat(ii).std=std(maxDat);
    sumdat(ii).min=min(maxDat);
    sumdat(ii).max=max(maxDat);
    sumdat(ii).median=median(maxDat);
    sumdat(ii).logmean=mean(log10(maxDat));
    sumdat(ii).logstd=std(log10(maxDat));
    
end

%% 

for ii=figList
    h=findobj(ii);
    h=findobj(ii,'type','axes');
    for jj=1:numel(h)
    boxY=h(jj).YLim;
    plot(h(jj),[1 1]*xLe,boxY,'LineStyle','--','Color',[0.3 0.3 0.3]);
    plot(h(jj),[1 1]*xTe,boxY,'LineStyle','--','Color',[0.3 0.3 0.3]);
    end
end

%%
for ii=figList
    h=findobj(ii);
    h=findobj(ii,'type','axes');
    boxY=h.YLim;
    boxX=h.XLim;
    t(1)=text(mean([xTe xLe]),min(boxY)/1.6,char(' Lower','Surface'));
    t(2)=text(mean([min(boxX) xLe]),min(boxY)/1.6,char('Upper','Surface'));
    %t(3)=text(mean([max(boxX) xTe]),min(boxY)/2.3,char(' Upper','Surface'));
    t(3)=text( xTe,min(boxY)/1.6,char('Trailing','   Edge'));
    t(4)=text( xLe,min(boxY)/1.6,char('Leading','  Edge'));
end
[t.HorizontalAlignment]=deal('center');
[t.Interpreter]=deal('latex');
[t.BackgroundColor]=deal([1 1 1]);
[t.FontSize]=deal(12);


%%


for ii=figList
    
     h=findobj(ii,'type','axes');
    if numel(h)>0
        for jj=1:numel(h)
           boxAxX=h(jj).XLim;
            boxAxY=h(jj).YLim;
            h(jj).NextPlot='add';
            for kk=[-19:2:19]
                l=plot(h(jj),[kk kk],[-1 2],'Color',[0.3 0.3 0.3],'Linestyle','--');
                
            end
            l.DisplayName='VOS Cell boundaries';
            
            c(ii).l=findobj(h(jj),'type','line');
            c(ii).l=c(ii).l(~cellfun(@isempty,{c(ii).l(:).DisplayName}));
            
            h(jj).XLim=boxAxX;
            h(jj).YLim=boxAxY;
        end
    end
end
%%

for ii=[1:5]
    h=findobj(ii,'type','axes');
    h.XLim=[-12 12];
end

%%

for ii=[2 6 7];
    
    h=findobj(ii,'type','axes');
    if numel(h)>0
        for jj=1:numel(h)
            lh=legend(flip(c(ii).l));
           lh.Interpreter='latex';
        end
    end
    
end

