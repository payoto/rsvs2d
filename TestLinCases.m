

rootDir='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Development Notes\TestLinOpt';
A=0.01:0.005:0.15;
for ii=1:numel(A)
    mkdir([rootDir,filesep,'test',int2str(ii)])
    [cdAll(ii),allobj(ii).obj]=SupersonicOptimLinRes(paramoptim,[rootDir,filesep,'test',int2str(ii)],-0.2,0.8,A(ii),500);
end
save(rootDir,'allobj','A')