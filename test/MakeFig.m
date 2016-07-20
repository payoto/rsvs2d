plotPoints= @(points) plot(points([1:end,1],1),points([1:end,1],2));
h=figure;
ax=axes;
hold on
ax.FontSize=12;
ax.TickLabelInterpreter='Latex';
l(1)=plotPoints(ZONE1);
l(2)=plotPoints(ZONE2);
for ii=1:2, l(ii).LineWidth=1;end
h.PaperPositionMode='auto';
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
legend(l(1:2),{'Optimised Profile','Parabolic Profile'},'interpreter','latex','location','SouthEast')
%%
pathStr='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Development Notes\Optimisation Results\CG\Interm\Dir_2016-05-15T010805_Test_CG_Aero\surf_fig.fig';
hgsave(1,pathStr);
pathStr2=pathStr;
pathStr2(end-2:end)='eps';
print(1,pathStr2,'-depsc');