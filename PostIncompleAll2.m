
n = 2;
for ii=1:n

    fid=fopen('PostDir2.txt','r');
    for jj=1:ii
        pathStr=fgetl(fid);
        
    end
    optstructPath = FindDir(pathStr, 'OptimRes', 0);
    nocalc=false;
    if numel (optstructPath)==1
        optstructPath = optstructPath{1};
    elseif numel(optstructPath)>0
        optstructPath2 = optstructPath(cellfun(@isempty, ...
            regexp(optstructPath,'_partial')));
        if numel(optstructPath)==0
            optstructPath2=optstructPath;
        end
        optstructPath = optstructPath2{1};
    else
        nocalc=true;
    end
    if ~nocalc 
        load(optstructPath)
        jj=0;
        while isempty(optimstruct(end-jj).population(1).objective)
            jj=jj+1;
        end
        optimstruct=optimstruct(1:end-jj-1);
         try
             disp(pathStr)
             disp(optstructPath)
            PostTreatIncomplete(pathStr,[],optimstruct);
         catch MEid
             T{ii}=MEid;

         end
    end
end
