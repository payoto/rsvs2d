% plotModes and deriv

%%
plotPoints= @(points) plot(points([1:end],1),points([1:end],2:end));
logPoints=@(points) plot(points(:,1),log10(abs(points(:,2:end))));

%%
kk=1;
pathStr{kk}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T115528_TestMode_sens_HF_0';kk=kk+1;
pathStr{kk}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T120450_TestMode_snak';kk=kk+1;
pathStr{kk}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T140306_TestMode_sens_HL_0';kk=kk+1;
pathStr{kk}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T142246_TestMode_sens_HL_1';kk=kk+1;
pathStr{kk}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T143106_TestMode_sens_HF_1';kk=kk+1;
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T152004_TestMode_snak_sens_HF_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T155616_TestMode_snak_snak_HFHA_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_07\Day_2017-07-31\Dir_2017-07-31T161204_TestMode_snak_sens_HFHA_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_08\Day_2017-08-01\Dir_2017-08-01T172417_TestMode_snak_sens_HF_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_08\Day_2017-08-01\Dir_2017-08-01T174502_TestMode_snak_snak_HF_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_08\Day_2017-08-02\Dir_2017-08-02T085803_TestMode_snak_sens_HFHA_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_08\Day_2017-08-02\Dir_2017-08-02T090818_TestMode_snak_snak_HFHA_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_08\Day_2017-08-02\Dir_2017-08-02T093846_TestMode_snak_sens_HF_';
pathStr{end+1}='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\Archive_2017_08\Day_2017-08-02\Dir_2017-08-02T094956_TestMode_snak_snak_HF_';

%%

for ii =1:numel(pathStr)
    [~,~,modecell{ii}]=ShowAreaError('foldermode',0,pathStr{ii},1,20);
end
[modecell{2,1:3}]=deal(14);
[modecell{2,4:end}]=deal(15);
for ii=1:size(modecell,2);
    modecell{3,ii}=modecell{1,ii}{2,modecell{2,ii}};
end

%%

close all
for ii=1:size(modecell,2);
    h(ii)=figure; subplot(121);
    plotPoints(FDPts(modecell{3,ii},3,1)),
    ax(2)=subplot(122);
    logPoints(FDPts(modecell{3,ii},3,1));
    ax(2).YLim=[-15 5];
    h(ii).Name=regexprep(pathStr{ii},'^.*TestMode_','');
end