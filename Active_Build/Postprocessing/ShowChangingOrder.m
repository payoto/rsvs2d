function h=ShowChangingOrder(ASOstruct,nameRun,beautify)
    
    if nargin<2
        nameRun='';
    end
    if nargin<3
        beautify=false;
    end
    
    if beautify
        for ii=1:numel(ASOstruct)
            ASOstruct(ii).obj=cummin(ASOstruct(ii).obj);
        end
    end
    
    fieldName='obj';
    nDat=cellfun(@numel,{ASOstruct.(fieldName)});
    boxDat=zeros([max(nDat),numel(ASOstruct)]);
    
    for ii=1:numel(ASOstruct)
        boxDat(1:nDat(ii),ii)=ASOstruct(ii).(fieldName);
        boxDat(nDat(ii)+1:end,ii)=boxDat(nDat(ii),ii);
    end
    ord=repmat(1:size(boxDat,2),[size(boxDat,1),1]);
    for ii=1:size(boxDat,1)
        [~,i]=sort(boxDat(ii,:));
        ord(ii,i)=ord(ii,:);
    end
    for ii=1:numel(ASOstruct)
        
        bestObj(ii)=(ASOstruct(ii).obj(end));
    end
    
    [~,i]=sort(ord(end,:));
    ord=ord(:,i);
    bestObj=flip(bestObj(i));
    
    ord=flip(ord,2);
    isLog=true;
    
    cDat=(bestObj);%ord(end,:);
    h=figure('Name',['Order Profile - ' nameRun]);
    
    if isLog
        cDat=log10(cDat);
    end
    cBounds=[min(cDat(isfinite(cDat))),max(cDat(isfinite(cDat)))];
    map=Viridis(2670);
    [datCol]=ProjectColormap(map,cDat,cBounds);
    l=plot(ord);
    [l.LineWidth]=deal(1);
    for ii=1:numel(l)
        l(ii).Color=datCol(ii,:);
    end
    
    xlabel('SNOPT major iteration')
    ylabel('Current Rank')
    
    if isLog
        h.Colormap=map(round(logspace(0,log10(size(map,1)-1),100)),:);
        h.Colormap=map(1:10:end,:);
    else
        h.Colormap=map(1:10:end,:);
    end
    ax=findobj(h,'type','axes');
    c=colorbar('peer',ax(1));
    
    if isLog
        c.Limits=cBounds;
        ax.CLim=c.Limits;
%         t=logspace(floor(cBounds(1)),ceil(cBounds(2)),1+10*(ceil(cBounds(2))-floor(cBounds(1))));
        t=[];
        for ii=floor(cBounds(1)):floor(cBounds(2))
            t=[t,(1:9)*10^ii];
        end
%         for ii=1:numel(t)
%             t(ii)=ceil(t(ii)*10^abs(floor(log10(t(ii)))))/10^abs(floor(log10(t(ii))));
%         end
        t=log10(unique(t));
        c.Ticks=t(t>=cBounds(1) & t<=cBounds(2));
        
        tickLabs=cell(size(c.Ticks));
        subTicks=[1 2 5];
        for ii=1:numel(c.Ticks);
            firstNum=round(round(10^c.Ticks(ii),1,'significant')/10^(floor(c.Ticks(ii))));
            
            if ii==1 || ii== numel(c.Ticks) ||  any(firstNum==subTicks)
                tickLabs{ii}=['$',int2str(firstNum),'\times10^{',int2str(floor(c.Ticks(ii))),'}$'];
            else
                tickLabs{ii}='';
            end
        end
        c.TickLabels=tickLabs;
        
        c.TickLabelInterpreter='latex';
    else
        c.Limits=cBounds;
        ax.CLim=c.Limits;
        c.TickLabelInterpreter='latex';
    end
    
    c.Label.String='Final $C_D$';
    c.Label.Interpreter='latex';
    c.Label.FontSize=ax(1).FontSize*ax(1).LabelFontSizeMultiplier;
    ax(1).XLim=[0,size(ord,1)];
    
end