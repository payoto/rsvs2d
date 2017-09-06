MoveToDir('source',1)
InitialiseSnakeFlow;
addpath test
ExecInclude


numList=[150,200,300,400,500];

for ii=1:numel(numList)
    [~,~,~,~,~,rootDir{ii}]=Main(sprintf('SnakNaca0012(%i)',numList(ii)));
    cellStr{ii}=num2str(numList(ii));
end
for ii=1:numel(numList)
    pathRestart=FindDir(rootDir{ii},'restart',0);
    load(pathRestart{1})
    grid.refined=gridrefined;grid.base=unstructReshape;grid.connec=connectstructinfo;snaxel=snakrestart.snaxel;
    grid2{ii}=grid;
    snax2{ii}=snaxel;
end
[nurbstruct,loop]=NURBSEngine('snakedense',snax2,grid2,cellStr);
h=findobj('type','figure');

QuickFigSave([h.Number]);