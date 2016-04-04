function [fill]=test_Convergence(snakSave)

%% Volume
errorFill=zeros([length(snakSave) length(snakSave(1).volumefraction)]);
for ii=1:length(snakSave),
    errorFill(ii,1:length(snakSave(ii).volumefraction))=([snakSave(ii).volumefraction(:).targetfill]-[snakSave(ii).volumefraction(:).volumefraction])./[snakSave(ii).volumefraction(:).targetfill];
end
fill=[snakSave(1).volumefraction(:).targetfill];
statError=[mean(errorFill);std(errorFill);mean(errorFill)-3*std(errorFill);mean(errorFill)+3*std(errorFill)];

figure
subplot(2,2,1)
plot(errorFill)
fftErr=fft(errorFill);
ylabel('fill Evolution')
subplot(2,2,2)
%plot(abs(fftErr))
plot(log10(abs(errorFill)))
ylabel('fill error')
%% Velocity

errorVel=zeros([length(snakSave) length(snakSave(1).volumefraction)]);
for ii=1:length(snakSave),
    
    newCellInd=[snakSave(ii).cellCentredGrid(:).index];
    meanSnaxVel=zeros(size(newCellInd));
    numSnax=zeros(size(newCellInd));
    for jj=1:length(newCellInd)
        meanSnaxVel(jj)=0;
        numSnax(jj)=numel(snakSave(ii).cellCentredGrid(jj).snaxel);
        if numSnax(jj)>0
            meanSnaxVel(jj)=sum(abs([snakSave(ii).cellCentredGrid(jj).snaxel(:).v]));
        
        end
    end
    cellNumSnax=zeros(size(snakSave(ii).volumefraction));
    cellSnaxVel=zeros(size(snakSave(ii).volumefraction));
    for jj=1:length(snakSave(ii).volumefraction)
        newSub=FindObjNum([],snakSave(ii).volumefraction(jj).newCellInd,newCellInd);
        cellSnaxVel(jj)=sum(meanSnaxVel(newSub));
        cellNumSnax(jj)=sum(numSnax(newSub));
    end
    cellNumSnax(cellNumSnax==0)=1e-6;
    errorVel(ii,1:length(snakSave(ii).volumefraction))=cellSnaxVel./cellNumSnax;
    
end

statVel=[mean(errorVel);std(errorVel);mean(errorVel)-3*std(errorVel);mean(errorVel)+3*std(errorVel)];

subplot(2,2,3)
plot(errorVel)
fftErr=fft(errorVel);

ylabel('Mean Absolute Velocity evolution')
subplot(2,2,4)
%plot(abs(fftErr))
plot(log10(abs(errorVel)))

ylabel('Logarithm of the mean absolute velocity')
end