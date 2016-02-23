function []=runSAVE(data,currentDir,saveDir,saveName,figureSave,struct)
% This function saves the current data and figure to a separate directory
% for later analysis.
% data: Is the variable to be saved
% currentDir: is the name of the current Matlab directory
% saveDir: is the substring defining the subdirectory in which the run
%           should be saved
% saveName: Is the file name to be given to the data
% struct: Is an optional argument to save the structre DATA as simple
%          variables in the save file

if ~exist('struct','var'); struct=0;end;
if ~exist('figureSave','var'); figureSave=0;end;

dir=[currentDir,saveDir];

system(['md "',dir,'"']);


if figureSave
    figh=1:20;
    fighindex=ishghandle(figh);
    figh=figh(fighindex);


    for ii=1:length(figh)

        saveas(figh(ii),[dir,'\figure',num2str(ii)],'fig');

    end
end

if struct
    
    save([dir,'\',saveName],'-struct','data')
    
else
   
    save([dir,'\',saveName],'data')
    
end


