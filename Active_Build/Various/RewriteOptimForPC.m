function [optimstruct]=RewriteOptimForPC(optimstruct,rootDir)
    
    
    for ii=1:numel(optimstruct)
        for jj=1:numel(optimstruct(ii).population)
            if ~isempty(optimstruct(ii).population(jj).location)
                [optimstruct(ii).population(jj)]=TrimProfLocation(optimstruct(ii).population(jj),rootDir);
            end
        end
    end
    
    
end



function [profile]=TrimProfLocation(profile,rootDir)
    
    
    posDir=regexp(profile.location,'Dir_');
    profile.location=[rootDir,filesep,profile.location(posDir:end)];
    profile.location=MakePathCompliant(profile.location);
    
    
end