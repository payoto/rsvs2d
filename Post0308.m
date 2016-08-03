

InitialiseSnakeFlow

kk=1;
pathStr{kk}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_07/Day_2016-07-31/Dir_2016-07-31T102545_bp3_Aero_CG_10_smooth_polypeak';
iterN(kk)=35;
kk=kk+1;
pathStr{kk}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_07/Day_2016-07-31/Dir_2016-07-31T135942_bp3_Aero_CG_10L_popk';
iterN(kk)=28;
kk=kk+1;
pathStr{kk}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_07/Day_2016-07-31/Dir_2016-07-31T145528_bp3_Aero_CG_10L_none';
iterN(kk)=40;
kk=kk+1;
pathStr{kk}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_07/Day_2016-07-31/Dir_2016-07-31T145559_bp3_Aero_CG_10L_pk';
iterN(kk)=40;
kk=kk+1;
pathStr{kk}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_07/Day_2016-07-31/Dir_2016-07-31T123053_bp3_Aero_CG_8_smooth_polypeak';
iterN(kk)=35;
kk=kk+1;

for ii=1:length(pathStr)
    PostTreatIncomplete(pathStr{ii},[1 iterN(ii)])
end