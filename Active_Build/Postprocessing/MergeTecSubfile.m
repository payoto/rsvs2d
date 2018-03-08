function []=MergeTecSubfile(tecSub,surfPlt)
    
    
    
    fidOrig=fopen(tecSub,'r');
    fidSurf=fopen(surfPlt,'r');
    if feof(fidOrig) || feof(fidSurf)
        error('Could not merge tecsubfiles')
    end
    % get headers
    for ii=1:9
        cellHeader{ii}=fgetl(fidOrig);
    end
    
    for ii=1:2
        cellStr{ii}=fgetl(fidSurf);
    end
    nElem=regexprep(regexp(cellStr{2},'E=[0-9]+','match'),'N=','');
    nElem=nElem{1};
    nNodes=regexprep(regexp(cellStr{2},'N=[0-9]+','match'),'N=','');
    nNodes=nNodes{1};
    
    
    cellHeader=regexprep(cellHeader,'STRANDID=3','STRANDID=5');
    cellHeader=regexprep(cellHeader,'NODES=[0-9]*',['NODES=',nNodes]);
    cellHeader=regexprep(cellHeader,'ELEMENTS=[0-9]*',['ELEMENTS=',nNodes]);
    
    
    nVarAdd=numel(regexp(cellHeader{1},','))-numel(regexp(cellStr{1},','));
    cellData=cell([1 str2num(nNodes)*2]);
    kk=1;
    catStr=[' ',int2str(zeros([1 nVarAdd]))];
    while ~feof(fidSurf)
        cellData{kk}=[fgetl(fidSurf),catStr];
        if kk==str2num(nNodes)
            catStr='';
        end
        kk=kk+1;
    end
    
    
    fclose(fidOrig);
    fclose(fidSurf);
    cellToAdd=[cellHeader,cellData];
    
    fidAdd=fopen(tecSub,'a');
    fprintf(fidAdd,'\n');
    WriteToFile(cellToAdd,fidAdd)
    
    fclose(fidAdd);
    % Define Correct Header
    % Modify strandId
    % modify number of elements
    
    % Expand data to be the same size
    % Concatenate files
    
    
    
    
    
    
end