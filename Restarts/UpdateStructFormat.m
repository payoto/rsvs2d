function iterstruct3=UpdateStructFormat(olditerstruct,newiterstruct)
    
    for ii=1:length(olditerstruct)
        
        for jj=1:length(olditerstruct(ii).population)
            fieldNames=fieldnames(olditerstruct(ii).population);
            for kk=1:length(fieldNames)
                
                newiterstruct(ii).population(jj).(fieldNames{kk})=olditerstruct(ii).population(jj).(fieldNames{kk});
                
            end
            
            newiterstruct(ii).population(jj).optimdat.var=1:length(newiterstruct(ii).population(jj).fill);
            newiterstruct(ii).population(jj).optimdat.value=newiterstruct(ii).population(jj).fill-newiterstruct(ii).population(1).fill;
        end
        
        
    end
    iterstruct3=newiterstruct;
    
    
end