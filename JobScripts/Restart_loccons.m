function Restart_loccons(ii,setN)
    MoveToDir('source',1)
    InitialiseSnakeFlow;
    
    
    restartDir='/panfs/panasas01/aero/ap1949/SnakVolParam/restarts/loccons/';
    if setN~=0
        restartDir=[restartDir,'/set',int2str(setN),'/'];
    end
    
    if ii~=0
        distinct=regexprep(num2str(ii),'\.','_');
        preStr='^LocConsTopo_prof_';
        postStr='__[0-9].mat';
        funcCall=['LocConsTopo_prof(',num2str(ii,'%e'),')'];
    else
       distinct='_Aero_DE_';
       preStr='bp3';
       postStr='_[0-9].mat';
       if setN==0
           funcCall=['bp3_Aero_DE_missile_horz'];
       else
           
           funcCall=['bp3_Aero_DE_smile_horz'];
       end
    end
    restartDir
    distinct
    preStr
    postStr
    [restartPath]=IdentifyRestart(restartDir,distinct,preStr,postStr);
    disp(funcCall)
    disp(restartPath)
    %ExecuteOptimisation(funcCall,{restartPath,{'DE',true}});
    

end