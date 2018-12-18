

pathToConstr='C:\Users\ap1949\Local Documents\PhD\res\aso_1804\aeroconv4';

pathDir=FindDir(pathToConstr,'Dir',1);

for ii=1:numel(pathDir);
    temp=FindDirRegex(pathDir{ii},'ASOper.*\.mat',0);
    perfPath{ii}=temp{1};
end

listInd=[1:21];
clear strctin; 
kk=1;
for ii=[listInd];
    strctin(kk)=load(perfPath{ii});
    kk=kk+1;
end
[~,c]=min(cellfun(@(x)numel(fieldnames(x)),{strctin.ASOstruct}));
fields=fieldnames(strctin(c).ASOstruct);
for ii=1:numel(strctin)
    for jj=1:numel(fields)
        [strctin2(ii).ASOstruct(1:numel(strctin(ii).ASOstruct)).(fields{jj})]=...
            deal(strctin(ii).ASOstruct.(fields{jj}));
    end
end
ASOconstr=[strctin2.ASOstruct];

save('fig/ASOMSconv4_new.mat','ASOconstr')
ASOconstr=ASOPerformanceAPI(ASOconstr,[],'splitCase','RunName','figList',[]);