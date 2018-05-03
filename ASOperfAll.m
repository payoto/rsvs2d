minII=1;
maxII=8;
fid=fopen('DirList.txt','r');
for ii=1:maxII
    pathStr{ii}=fgetl(fid);
end
fclose(fid);

for ii=minII:maxII
    ASOstructCell{ii}=ASOPerformanceAPI(pathStr{ii},2);
end

ASOstruct2=[ASOstructCell{:}];
ASOPerformanceAPI(ASOstruct2,1,'dirSave','fig','nameRun','lvl 1-5 errmode all')


