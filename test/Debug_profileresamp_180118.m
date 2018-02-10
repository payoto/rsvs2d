%%
InitialiseSnakeFlow
ExecInclude
InverseDesign_Error
%%
global kk
kk=1;
iter50_2=load('C:\Users\ap1949\Local Documents\PhD\res\tri_180111\Dir_2018-01-16T051943_volsweeplocal_9_000000e-02_uov_tri_1_2_3\profile_50_2\restart_50_2.mat')
iter51_1=load('C:\Users\ap1949\Local Documents\PhD\res\tri_180111\Dir_2018-01-16T051943_volsweeplocal_9_000000e-02_uov_tri_1_2_3\profile_51_1\restart_51_1.mat')
iter51_3=load('C:\Users\ap1949\Local Documents\PhD\res\tri_180111\Dir_2018-01-16T051943_volsweeplocal_9_000000e-02_uov_tri_1_2_3\profile_51_3\restart_51_3.mat')
iter(1)=iter51_1;
iter(2)=iter50_2;
iter(3)=iter51_3;
%%
for ii=1:numel(iter)
    iter(ii).loop.snaxelcoord=iter(ii).loop.snaxel.coord;
end
paramspline.splineCase='smoothpts';
for ii=1:numel(iter);
    iter(ii).loop.subdispline2=ResampleSpline(iter(ii).loop.subdivision,paramspline);
end
%%
%close all
fields={'snaxelcoord','subdivision','subdivspline','subdispline2'};
for ii=1:numel(fields)
    analysisCoord=iter(1).loop.(fields{ii});
    for jj=2:3
        targCoord=iter(jj).loop.(fields{ii});
        analysisCoord=RemoveIdenticalConsecutivePoints(analysisCoord);
        targCoord=RemoveIdenticalConsecutivePoints(targCoord);
        [errorMeasure,modifiedDistance]=CompareProfilesDistance2(analysisCoord,targCoord);
        plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
        h=figure;
        subplot(2,1,1)
        plotPoints(analysisCoord)
        hold on
        plotPoints(targCoord)
        legend('snake points','Target points')
        ax=subplot(2,1,2);
        plotPoints(modifiedDistance)
        ax.YScale='log';
    end
end

%%
clear parList*
clear newPts*
clear normPoints*
clear newParList*
clear curvParam*
clear curvParamErr*
kk=1;
for ii=1:3;
for jj=1:numel(iter(ii).loop)
    iter(ii).loop(jj).subdispline2=ResampleSpline(iter(ii).loop(jj).subdivision,paramspline);
end
end
kk=4;
funcCell={@(v) rand(v) , @(v) ones(v) , @(v) zeros(v)};
for ii=1:3;
    iter2(ii).loop.subdispline2=ResampleSpline(iter(1).loop.subdivision+(1e-16+funcCell{ii}(size(iter(1).loop.subdivision))*2e-16),paramspline);
end
%%
figure, plot(parList1-parList4)
hold on
plot(parList1-parList5)
plot(parList1-parList6)
figure, hold on
plot(curvParam1)
plot(parList1)
figure, plot(curvParam1-curvParam4)
hold on
plot(curvParam1-curvParam5)
plot(curvParam1-curvParam6)
figure, plot(log10(abs(curvParam1-curvParam4)))
hold on
plot(log10(abs(curvParam1-curvParam5)))
plot(log10(abs(curvParam1-curvParam6)))
curvDiff=log10(abs(curvParam1-curvParam4));
curvDiffAve=zeros(size(curvDiff));
for ii=1:size(curvDiff,1);
    curvDiffAve(ii,:)=max(curvDiff(mod(ii-2:ii+1,size(curvDiff,1))+1,:));
end
curvDiffAve=max(curvDiffAve,[],2);
curvParamNorm1=sqrt(sum(curvParam1.^2,2));
plot(curvDiffAve)
figure,plot(-log10(abs(curvParamNorm1))+curvDiffAve)
figure,plot(log10(abs(curvParamNorm1)))
hold on
plot(curvDiffAve)
figure, plot(curvParamErr1)
hold on
plot(curvParamErr4)
plot(curvParamErr5)
plot(curvParamErr6)
figure, plot(log10(abs(curvParamErr1)))
hold on
plot(log10(abs(curvParamErr4)))
plot(log10(abs(curvParamErr5)))
plot(log10(abs(curvParamErr6)))

%% testing multi-body

ii=6
figure,
plotPoints= @(points,f) plot(points([1:end,1],1),points([1:end,1],2),['-',f]);
for jj=1:numel(iter(ii).loop)
    iter(ii).loop(jj).subdispline2=ResampleSpline(iter(ii).loop(jj).coord,paramspline);
end
figure, hold on
for jj=1:numel(iter(ii).loop)
    plotPoints(iter(ii).loop(jj).subdispline2,'o');
end
for jj=1:numel(iter(ii).loop)
    plotPoints(iter(ii).loop(jj).coord,'');
end