function [out]=SnakeValid(validationName)
    
    include_PostProcessing
    
    snakCell={'Snakestestsmooth1','Snakestestsmooth1_2','Snakestestsmooth2',...
        'Snakestestsmooth3','Snakestestsmooth3_1','Donught','Donught2',...
        'SnakesFoilVVSmall','BuzmanBiplane3','SnakesFoilVVSmall4',...
        'WeirdShapeIn','WeirdShapeOut'};
    
    T{length(snakCell)}=[];
    
    out=repmat(struct('unstructured',[],'loop',[],'unstructReshape',[],'snakSave',[]...
    ),[1 length(snakCell)]);
    procdat=repmat(struct('warningnum',0,'errornum',0,'cputime',0,'termination',false, 'niter',0,...
        'snaketime',0,'volerror',1,'velerror',1,'length',0),[1 length(snakCell)]);
    
    parfor ii=1:length(snakCell)
        pause(ii/5)
        [T{ii},out(ii)]=CallMain(snakCell{ii});
        [procdat(ii)]=ProcessData(out(ii),T{ii},procdat(ii));
    end
    
    
    
    
    
    
end

function [T,out]=CallMain(snakCell)
    
    [T,out.unstructured,out.loop,out.unstructReshape,out.snakSave]=...
    evalc('Main([''val_'',snakCell])');
end

function [procdat]=ProcessData(out,T,procdat)
    % out is the output structure containing the data 
    % T is the char array containing the text output
    
    
    
    try
        % read from text
        procdat.warningnum=numel(regexpi(T,'warning'));

        procdat.errornum=numel(regexpi(T,'error'));

        timeMark=regexp(regexpi(T,'Iteration time[^\n]*','match'),':.*','match');
        procdat.cputime=str2num(regexprep(timeMark{1}{1},':',' '))*[3600;60;1;0.001]; % time between start and end. (look for "Iteration Time:")

        procdat.termination=~isempty(regexpi(T,'Snakes Converged!')); % true / false (true means terminate with convergence) (look for Snakes Converged!)

        % read from array
        procdat.niter=length(out.snakSave); % number of iteration to termination

        procdat.snaketime=sum([out.snakSave(:).dt]);
        procdat.volerror=out.snakSave(end).currentConvVolume;
        procdat.velerror=out.snakSave(end).currentConvVelocity;
        procdat.length=out.snakSave(end).lSnak;
    catch
        procdat.errornum=1;
    end
    
end