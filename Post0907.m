ii=1;
pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-06/Dir_2016-09-06T142120_bp3_missile2';
ii=ii+1;
pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-06/Dir_2016-09-06T143402_bp3_MT_pop_25';
ii=ii+1;
pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-06/Dir_2016-09-06T143424_bp3_MT_pop_50';
ii=ii+1;
pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-06/Dir_2016-09-06T143449_bp3_MT_pop_75';
ii=ii+1;
pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-06/Dir_2016-09-06T143514_bp3_MT_pop_100';
ii=ii+1;
pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-06/Dir_2016-09-06T143612_bp3_MT_pop_125';
ii=ii+1;
pathStr{ii}='/panfs/panasas01/aero/ap1949/SnakVolParam/results/Optimisation/Archive_2016_09/Day_2016-09-06/Dir_2016-09-06T143640_bp3_MT_pop_150';

numIter=[12 48 28 19 14 11 9];
for ii=1:length(pathStr)
    PostTreatIncomplete(pathStr{ii},[1,numIter(ii)])
end