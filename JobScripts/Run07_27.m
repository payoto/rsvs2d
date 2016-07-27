MoveToDir('source',1)
InitialiseSnakeFlow



optimstruct=ExecuteOptimisation('desk_Aero_CG_20L_po');
save('optimPolyRes.mat','optimstruct')
ExecuteOptimisation('desk_Aero_CG_20L_pk',{'optimPolyRes',{'conjgrad'},true})

    



