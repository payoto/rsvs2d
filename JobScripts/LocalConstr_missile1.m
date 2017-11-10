
MoveToDir('source',1)
InitialiseSnakeFlow;

% try 
    ExecuteOptimisation('bp3_Aero_CG_missile_in');
% catch ME
%     ME.getReport
% end
try 
    ExecuteOptimisation('bp3_Aero_CG_missile_out');
catch ME
    ME.getReport
end