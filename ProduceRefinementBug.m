%% Load Files
InitialiseSnakeFlow
include_Mex_Wrapper

load('refinementbug.mat')


%% Truncate



%% Execute

[gridrefined,connecstructinfo]=GridRefine_Wrapper(gridIn,input2, input3);
