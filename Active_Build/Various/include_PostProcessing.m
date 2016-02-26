function [] = include_PostProcessing()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end

%%

%% General
function []=savefig(figh,figname)

print(figh,figname,'-djpeg','-r600')

end

function [marker,t]=GenerateResultMarker(typDat)
    
    t=now;
    marker=[datestr(t,'yyyy-mm-ddTHHMMSS')...
        ,'_',typDat];
    
end

function [resultDirectory]=GenerateResultDirectoryName(marker,resultRoot,archiveName,t)
    if ~exist('t','var'),t=now;end
    dateSubFolders=['Archive_',datestr(now,'yyyy_mm'),'\Day_',datestr(t,29)];
    resultDirectory=[resultRoot,'\',archiveName,'\',dateSubFolders,...
        '\','Dir_',marker];
    system(['md "',resultDirectory,'"']);
end

function []=WriteToFile(cellLoops,FID)
    % writes the data to the file
    
    
    for ii=1:length(cellLoops)
        if numel(cellLoops{ii})>0
            for jj=1:length(cellLoops{ii}(:,1))
                fprintf(FID,'%s \n',cellLoops{ii}(jj,:));
            end
        else
            fprintf(FID,'\n');
        end
    end
end

function []=ReorganiseSubPlots(ax,sizSubPlot,outerPad,interPad,fontSizes)
    
    if ~exist('outerPad','var'),outerPad=[0.1,0.1,0.05,0.05];end
    if ~exist('interPad','var'),interPad=[0.05,0.05];end
    if ~exist('fontSizes','var'),fontSizes=[10,12];end
    
    xHeight=1-(outerPad(1)+outerPad(2));
    xLength=(xHeight-interPad(1)*(sizSubPlot(1)-1))/sizSubPlot(1);
    xStart=outerPad(1);
    
    yHeight=1-(outerPad(3)+outerPad(4));
    yLength=(yHeight-interPad(2)*(sizSubPlot(2)-1))/sizSubPlot(2);
    yStart=1-outerPad(4)-yLength;
    
    for ii=1:length(ax)
        jj=ii;
        xInd=rem(jj-1,sizSubPlot(1));
        yInd=ceil(jj/sizSubPlot(1))-1;
        
        loXCorn=(xInd)*(xLength+interPad(1))+xStart;
        loYCorn=-(yInd)*(yLength+interPad(2))+yStart;
        
        posVec=[loXCorn,loYCorn,xLength,yLength];
        set(ax(ii),'position',posVec)
        set(ax(ii),'fontsize',fontSizes(1))
        set(get(ax(ii),'xlabel'),'fontsize',fontSizes(2))
        set(get(ax(ii),'ylabel'),'fontsize',fontSizes(2))
        set(get(ax(ii),'zlabel'),'fontsize',fontSizes(2))
        
    end
    
end

%% Parameter File

function []=GenerateParameterFile(FID,param,t,marker)
    
    paramCell{1}='% Parameter File';
    paramCell{2}=['% ',datestr(t)];
    paramCell{3}=['% ',marker];
    for ii=length(param.structdat.vars):-1:1
        [paramCell{ii+4}]=ExtractVariablePathAndValue(param,ii);
        
    end
    
    WriteToFile(paramCell,FID);
    fclose(FID);
end

function [paramStr]=ExtractVariablePathAndValue(param,varNum)
    
    varstruct=param.structdat.vars(varNum);
    actstruct=param;
    pathVar='param';
    for jj=1:length(varstruct.vec)
        pathVar=[pathVar,'.',param.structdat.fields{varstruct.vec(jj)}];
        actstruct=actstruct.(param.structdat.fields{varstruct.vec(jj)});
    end
    
    varVal=actstruct;
    [varStr]=GenerateVariableString(varVal);
    paramStr=[pathVar,' = ',varStr,';'];
end

function [varStr]=GenerateVariableString(startVar)
    
    classVar=class(startVar);
    [m,n]=size(startVar);
    switch classVar
        case 'char'
            
            openStr='[';
            closeStr=']';
            for ii=1:m
                varStrCell{ii,1}=['''',startVar(ii,:),''''];
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,1);
        case 'cell'
            
            openStr='{';
            closeStr='}';
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=GenerateVariableString(startVar{ii,jj});
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
        case 'double'
            
            openStr='[';
            closeStr=']';
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=num2str(startVar(ii,jj));
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
        case 'logical'
            openStr='[';
            closeStr=']';
            for ii=1:m
                for jj=1:n
                    if startVar(ii,jj) 
                        curStr='true';
                    else
                        curStr='false';
                    end
                    varStrCell{ii,jj}=curStr;
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
        otherwise
            if ~isempty(regexp(classVer,'int','once'))
                
            openStr='[';
            closeStr=']';
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=int2str(startVar(ii,jj));
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
            end
            warning('Class is not catered for and will not be printed correctly to parameter file')
    end
    
end

function [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n)

    
    varStr=openStr;
    for ii=1:m-1
        for jj=1:n-1
            varStr=[varStr,varStrCell{ii,jj},','];
        end
        varStr=[varStr,varStrCell{ii,n},';'];
    end
    for jj=1:n-1
        varStr=[varStr,varStrCell{m,jj},','];
    end
    varStr=[varStr,varStrCell{m,n},closeStr];
end

