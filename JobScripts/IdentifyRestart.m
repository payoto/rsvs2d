function [restartPath]=IdentifyRestart(restartDir,distinct,preStr,postStr)
    
    c=dir(restartDir);
    
    allNames={c(:).name};
    restartNamesFull=allNames(~cellfun(@isempty,regexp(allNames,'OptimRes')));
    
    restartNames=regexprep(restartNamesFull,preStr,'');
    restartNames=regexprep(restartNames,postStr,'');
    
    if isnumeric(distinct)
        restartNames=regexprep(restartNames,'_','.');
        restartNum=cellfun(@str2num,restartNames,'UniformOutput',false);
        
        ismatch=false;
        kk=0;
        while ~ismatch && kk<numel(restartNum)
            kk=kk+1;
            ismatch=all(distinct==restartNum{kk});
            if isempty(ismatch)
                ismatch=false;
            end
        end
        
        if ~ismatch 
            error('Restart not found')
        end
        
        restartPath=[restartDir,filesep,restartNamesFull{kk}];
        
    elseif ischar(distinct)
        
        kk=find(~cellfun(@isempty,regexp(restartNames,distinct)));
        
        
        if isempty(kk)
            error('Restart not found')
        end
        
        restartPath=[restartDir,filesep,restartNamesFull{kk(1)}];
    else
        error('Not supported')
    end
    
end