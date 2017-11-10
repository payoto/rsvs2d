function paramoptim=AddFieldToParam(fieldName,fieldVal,paramoptim)
    
    
    paramoptim.general.(fieldName)=fieldVal;
    nMatch=numel(paramoptim.structdat.vars)+1;
    paramoptim.structdat.vars(nMatch).name = fieldName;
    nField=numel(paramoptim.structdat.fields)+1;
    paramoptim.structdat.fields{nField}=fieldName;
    paramoptim.structdat.vars(nMatch).vec = [1 nField];
    paramoptim.structdat.vardat.varmatch(regexp(paramoptim.structdat.vardat.names,fieldName))=nMatch;
    paramoptim.structdat.vardat.names=[paramoptim.structdat.vardat.names,fieldName,'#'];
    paramoptim.structdat.vardat.varmatch(regexp(paramoptim.structdat.vardat.names,'rankType'))=nMatch;
    
    
end