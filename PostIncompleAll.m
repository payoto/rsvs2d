for ii=1:5

    fid=fopen('PostDir2.txt','r');
    for jj=1:ii
        pathStr=fgetl(fid);
        optstructPath=fgetl(fid);
    end
    load(optstructPath)
    
    jj=0;
    while isempty(optimstruct(end-jj).population(1).objective)
        jj=jj+1;
    end
    optimstruct=optimstruct(1:end-jj-1);
%     try
        PostTreatIncomplete(pathStr,[],optimstruct);
%     catch MEid
%         T{ii}=MEid;
%         
%     end
end
