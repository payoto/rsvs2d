minII=1;
maxII=6;
fid=fopen('DirList.txt','r');
pathStr=cell([maxII,2]);
for ii=1:maxII
    pathStr(ii,1:2)=regexp(fgetl(fid),'\s+','split');
end
fclose(fid);
T=cell([maxII,1]);
parfor ii=minII:maxII
    try 
    ASOPerformanceAPI(pathStr{ii,1},str2double(pathStr{ii,2}));
    catch MEid
        T{ii}=MEid;
    end
end



