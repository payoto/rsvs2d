%% Execute Mesh Generation

[unstructured,loop,unstructReshape,snakSave]=Main('WedgPlot')

% %% Execute grid refinement
% 
% [gridrefined,connectinfo]=GridRefinement(unstructReshape,8,true);
% %% Executes Snakes
% 
% [snaxel,snakposition,snakSave,loopsnaxel]=Snakes(unstructured,loop,5,0,50);
% 
