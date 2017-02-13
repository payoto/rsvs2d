% Script to introduce desirable modifications to ExecuteOptimisation

flagOut=false;
load('Debug_1_VELMIN')
baseGrid=grid.base
gridrefined=grid.refined
ncell=numel(grid.base.cell);
[connectstructinfo.cell(1:ncell).new]=deal(grid.connec.cell(:).newCellInd);
[connectstructinfo.cell(1:ncell).old]=deal(grid.connec.cell(:).oldCellInd);
connectstructinfo.edge=struct('new',[],'old',[]);
iterstruct=optimstruct;
maxIter=0;