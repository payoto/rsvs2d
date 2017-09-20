% small script to prepare data for saving


load GlobalRunIndex

NewRun=max(GlobalRunNum)+1;

GlobalRunNum=[GlobalRunNum,NewRun];

save(GlobalRunIndex,'GlobalRunNum')