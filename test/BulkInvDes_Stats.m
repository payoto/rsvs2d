function [dat]=BulkInvDes_Stats(optimstruct)
    
    dat=struct(...
        'Mean',[],...
        'reflevel',[],...
        'Median',[],...
        'max',[],...
        'satisfied',[],...
        'desVar',[],...
        'active',[]...
        );
    for ii=1:numel(optimstruct)
        obj=[optimstruct(ii).population(:).objective];
        dat(ii).Mean = mean(obj);
        dat(ii).reflevel = ii-1;
        dat(ii).Median = median(obj);
        dat(ii).max = max(obj);
        dat(ii).satisfied = sum(obj<8e-4)/numel(obj);
        
        pathStr=optimstruct(ii).population(1).location;
        posI = regexp(pathStr,'iteration');
        pathStr = pathStr(1:posI-1);
        
        currPopPath=FindDir(pathStr,'OptimRes',0);
        
        optstructCurr = load(currPopPath{1});
        
        dat(ii).desVar = numel(optstructCurr.optimstruct(ii).population(1).fill);
        dat(ii).active = 0;
        for jj = 1:numel(optimstruct(ii).population)
            dat(ii).active = dat(ii).active + ...
                sum(optstructCurr.optimstruct(ii).population(jj).fill~=0 ...
                    & optstructCurr.optimstruct(ii).population(jj).fill~=1);
        end
        dat(ii).active = dat(ii).active/numel(optimstruct(ii).population);
    end

end