load('C:\Users\ap1949\Local Documents\PhD\res\invdes_170403\3_uu_contcurvevol_0_1\Dir_2017-03-31T212022_3_uu_contcurvevol_0_1_3\iteration_121\profile_1\restart_121_1.mat')
restartsnake=restartsnak


load('C:\Users\ap1949\Local Documents\PhD\res\invdes_170403\3_uu_contcurvevol_0_1\Dir_2017-03-31T212022_3_uu_contcurvevol_0_1_3\OptimRes_2017-03-31T212022_invdeslocal_test3_uu_contcurvevol_0_1__1_2_3partial.mat')
[iterstruct(121).population(:).fill]=deal(optimstruct(121).population(:).fill);

paramoptim.parametrisation.snakes.step.maxStep=10;