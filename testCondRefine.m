oldGrid=SelectRefinementCells(optimstruct(end).population,grid,paramoptim)


refStep=1

varNames={'refineOptim'};
refineCellLvl=ExtractVariables(varNames,paramoptim);

refparamsnake=SetVariables({'refineGrid','typeRefine'},{refineCellLvl(refStep,:),'automatic'},...
paramoptim.parametrisation);

[~,baseGrid,gridrefined,connectstructinfo,~,~]=GridInitAndRefine(refparamsnake,oldGrid.base);

include_GridCheck

CheckGrid(gridrefined)