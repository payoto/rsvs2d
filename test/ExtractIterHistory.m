
function [dataTable]=ExtractIterHistory(rawOut)
    
    dat=regexp(rawOut,'\n','split');
    logArray=~cellfun(@isempty,regexp(dat,'^\s*[0-9][0-9e+\-\.\s]*$'));
    logMaxRes=~cellfun(@isempty,regexp(dat,'log10\[Maximum residual\]:'));
    datNum=str2num(char(dat{logArray})); %#ok<ST2NM>
    fields=regexp(dat{find(logArray,1,'first')-1},'\s*','split');
    fields=fields(~cellfun(@isempty,fields));
    
    for ii=1:size(datNum,2)
        fields{1,ii}=matlab.lang.makeValidName(fields{1,ii});
        fields{2,ii}=datNum(:,ii);
    end
    fields{1,ii+1}='maxResidual';
    fields{2,ii+1}=str2num(char(regexprep(regexprep(...
        dat(logMaxRes),'\.\s*$',''),'log10\[Maximum residual\]:',''))); %#ok<ST2NM>
    fields{1,ii+2}='maxResidualPos';
    fields{2,ii+2}=str2num(char(regexprep(regexprep(dat(find(logMaxRes)+1),...
        '^.*(',''),').*$',''))); %#ok<ST2NM>
    dataTable=struct(fields{:});
    
    
end