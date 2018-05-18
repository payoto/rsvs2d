function h=ShowChangingOrder(ASOstruct,nameRun)
    
    if nargin<2
        nameRun='Order Profile';
    end
    for jj=1:3
        for ii=1:numel(ASOstruct)
            ASOstruct(ii).obj([false;ASOstruct(ii).obj(2:end)>ASOstruct(ii).obj(1:end-1)])=...
                ASOstruct(ii).obj(find([false;ASOstruct(ii).obj(2:end)>ASOstruct(ii).obj(1:end-1)])-1);
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
    
    [~,i]=sort(ord(end,:));
    ord=ord(:,i);
    
    h=figure('Name',nameRun);
    [datCol]=ProjectColormap(h.Colormap,ord(end,:),[1,size(ord,2)]);
    l=plot(ord);
    
    for ii=1:numel(l)
        l(ii).Color=datCol(ii,:);
    end
        
end