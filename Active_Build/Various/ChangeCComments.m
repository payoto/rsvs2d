function []=ChangeCComments()
    
    filePath='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\MEX_Function_Directory\MEX_Sources\findobjnum';
    fileList={'FindObjNum_MEX.c'};
    
    for ii=1:length(fileList)
        RewriteFile(filePath,fileList{ii})
    end
    
end


function []=RewriteFile(filePath,fileName)
    
    fIDR=fopen([filePath,filesep,fileName],'r');
    fIDW=fopen([filePath,filesep,'Modif_',fileName],'w');
    
    isInComment=0;
    while ~feof(fIDR)
        lineStr=fgetl(fIDR);
        lineStr=regexprep(lineStr,'%','%%');
        lineStr=regexprep(lineStr,'\\','\\\\');
        comStart=regexp(lineStr,'//');
        isInComment=isInComment+(~isempty(regexp(lineStr,'/\*', 'once')))...
            -(~isempty(regexp(lineStr,'*/', 'once')));
        if ~isempty(comStart)
            lineStr(comStart:comStart+1)='/*';
            lineStr(end+1:end+3)=' */';
            if isInComment
                lineStr(end+1:end+3)=' /*';
            end
        end
        
        fprintf(fIDW,[lineStr,'\n']);
    end
    
    fclose(fIDR);
    fclose(fIDW);
    
end


