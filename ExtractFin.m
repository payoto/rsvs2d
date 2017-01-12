pathStr={'./Dir_2016-12-23T165538_areabusesweep_3_000000e-02_/OptimRes_2016-12-23T165538_areabusesweep_3_000000e-02_.mat'
'./Dir_2016-12-23T165601_areabusesweep_2_000000e-02_/OptimRes_2016-12-23T165601_areabusesweep_2_000000e-02_.mat'
'./Dir_2016-12-23T165608_areabusesweep_4_000000e-02_/OptimRes_2016-12-23T165608_areabusesweep_4_000000e-02_.mat'
'./Dir_2016-12-24T090549_areabusesweep_5_000000e-02_/OptimRes_2016-12-24T090549_areabusesweep_5_000000e-02_.mat'
'./Dir_2016-12-24T090622_areabusesweep_6_000000e-02_/OptimRes_2016-12-24T090622_areabusesweep_6_000000e-02_.mat'
'./Dir_2016-12-24T093728_areabusesweep_7_000000e-02_/OptimRes_2016-12-24T093728_areabusesweep_7_000000e-02_.mat'
'./Dir_2016-12-24T093805_areabusesweep_8_000000e-02_/OptimRes_2016-12-24T093805_areabusesweep_8_000000e-02_.mat'
'./Dir_2016-12-24T130246_areabusesweep_9_000000e-02_/OptimRes_2016-12-24T130246_areabusesweep_9_000000e-02_.mat'
'./Dir_2016-12-24T130326_areabusesweep_1_000000e-01_/OptimRes_2016-12-24T130326_areabusesweep_1_000000e-01_.mat'
'./Dir_2016-12-24T143129_areabusesweep_1_100000e-01_/OptimRes_2016-12-24T143129_areabusesweep_1_100000e-01_.mat'
'./Dir_2016-12-24T143207_areabusesweep_1_200000e-01_/OptimRes_2016-12-24T143207_areabusesweep_1_200000e-01_.mat'
'./Dir_2017-01-01T163609_areabusesweep_5_000000e-02_/OptimRes_2017-01-01T163609_areabusesweep_5_000000e-02_.mat'
'./Dir_2017-01-01T163614_areabusesweep_7_000000e-02_/OptimRes_2017-01-01T163614_areabusesweep_7_000000e-02_.mat'
'./Dir_2017-01-01T163640_areabusesweep_6_000000e-02_/OptimRes_2017-01-01T163640_areabusesweep_6_000000e-02_.mat'
'./Dir_2017-01-01T163650_areabusesweep_9_000000e-02_/OptimRes_2017-01-01T163650_areabusesweep_9_000000e-02_.mat'
'./Dir_2017-01-01T163655_areabusesweep_1_100000e-01_/OptimRes_2017-01-01T163655_areabusesweep_1_100000e-01_.mat'
'./Dir_2017-01-01T163726_areabusesweep_1_000000e-01_/OptimRes_2017-01-01T163726_areabusesweep_1_000000e-01_.mat'
'./Dir_2017-01-01T163733_areabusesweep_1_200000e-01_/OptimRes_2017-01-01T163733_areabusesweep_1_200000e-01_.mat'};
path2='C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\results\Optimisation\HPC\DE_population\'
for ii=1:numel(pathStr)
    load([path2,pathStr{ii}])
    [minVal(ii),i]=min([optimstruct(end).population(:).objective]);
    minA(ii)=optimstruct(end).population(i).additional.A;
end