

paramoptim.parametrisation = SetVariables({'modeSmoothType'},{'optimlinprog'},paramoptim.parametrisation);

paramsnake = paramoptim.parametrisation;
[out]=SnakesSensitivity('fill',grid.refined,restartsnak,grid.base,paramoptim.parametrisation,paramoptim,rootFill, out)
[out]=SnakesSensitivity('profile',grid.refined,restartsnak,grid.base,paramoptim.parametrisation,paramoptim,rootFill, out)

figure, hold on,
for ii = 1:numel(out)
    x = [1:numel(out(1).loopsens.snaxel.coord(:,1))]';
    y = zeros(size(out(1).loopsens.snaxel.coord(:,1)));
    u = out(ii).loopsens.snaxel.coord(:,1)-out(1).loopsens.snaxel.coord(:,1);
    v = out(ii).loopsens.snaxel.coord(:,2)-out(1).loopsens.snaxel.coord(:,2);
    dotprod=sum(normDirection.*[u,v],2);
    [~,imax] = max(abs(dotprod));
    plot(x,((dotprod/dotprod(imax))))
end