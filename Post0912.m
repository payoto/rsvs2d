ii=1;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T022451_bp3_missile2';
ii=ii+1 ;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T032754_bp3_smile';
ii=ii+1 ;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T032846_bp3_MT_pop_25';
ii=ii+1 ;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T094131_bp3_MT_pop_75';
ii=ii+1 ;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T094151_bp3_MT_pop_100';
ii=ii+1 ;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T100308_bp3_MT_pop_125';
ii=ii+1 ;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T101128_bp3_MT_pop_50';
ii=ii+1 ;
 pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-07/Dir_2016-09-07T101309_bp3_MT_pop_150';

numIter=[36 34 142 50 37 28 70 23];
for ii=1:length(pathStr)
    PostTreatIncomplete(pathStr{ii},[1,numIter(ii)])
end