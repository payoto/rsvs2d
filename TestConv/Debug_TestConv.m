
InitialiseSnakeFlow
ExecInclude

addpath('/panfs/panasas01/aero/ap1949/ASO_LK/SRC')
addpath('/panfs/panasas01/aero/ap1949/ASO_LK/SRC/matlab-snopt')
 
 profVecNoAdj=[62,72,65,66,77,51,8]; % all sensible
 profVecNoDir=[37,20,33,40,50]; % 33-50 are weird
 profVecWorks=[4 10 13 16 ]; % all sensible
 pathStr='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2018_04/Day_2018-04-23/Dir_2018-04-23T102550_ASOMS_subdiv_5_basis_1_';
 cflNum=[1 2 4];
 for ii=1:numel(cflNum)
    [ASOConvstruct,population]=RunASOFlowConvTest(pathStr,...
        '/panfs/panasas01/aero/ap1949/SnakVolParam/results/ConvSU2',2,...
        [profVecNoAdj,profVecNoDir,profVecWorks],['CFL(',num2str(cflNum(ii)),')']);
 end