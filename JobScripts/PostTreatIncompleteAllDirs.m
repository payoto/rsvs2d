
MoveToDir('source',1)
InitialiseSnakeFlow;


 [T,postDir]=IdentifyPostTreatTargets({'/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation'},4,[]);

for ii=1:numel(T)
	
	fprintf('\n')
	if ischar(T{ii})
		fprintf(postDir{ii})
		fprintf(T{ii})
		fprintf('\n')
	else
		fprintf(postDir{ii})
		
		T{ii}.getReport
	end
end
