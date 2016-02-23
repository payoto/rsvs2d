%% Post Treatment of result All in 1
% This script is used to find different runs in a directory extract 
% information about them and compute the RMSE of the surrogate models
% against a base model.

%% Workspace preparation
close all
clear all
clc

%% Directory Load and Folder Operations

% Root folder to be analysed

root_dir='C:\Users\Alexandre\Documents\Travail\Uni\IXP\Simulations\XFOIL6.99\3D_200points\';

d = dir(root_dir);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
logiFolds=regexp(nameFolds,'run');



for ii=1:length(nameFolds)
    if isempty(logiFolds{ii})
        logiFolds{ii}=0;
    end
end
logiFolds=logical(cell2mat(logiFolds));

nameFolds=nameFolds(logiFolds);
% numFolds=regexprep(nameFolds,'run','');
% numFolds=cell2mat(numFolds);


n_run=sum(logiFolds);

%% Load reference data

load([root_dir,'bigdata.mat'])

data_ref=data;

clear data

%% Analysis of Subfolders

adap_var_array=[1, 2 ,3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3];

for ii=1:n_run
    
    %% Subdirectory specific variable extraction
    
    fprintf('Run %i out of %i started\n',ii,n_run)
    
    run_dir=[root_dir,nameFolds{ii},'\'];
   
    load([run_dir,'RunParam.mat'])
    
    run_param=data;
    clear data
    
    post(ii).desVar=run_param.desVar;
    desVar=run_param.desVar;
    post(ii).des_constant=run_param.des_constant;
    des_constant=run_param.des_constant;
    post(ii).N_LHS=run_param.N_LHS;
    N_LHS=run_param.N_LHS;
    post(ii).N_adapCycles=run_param.N_adapCycles;
    N_adapCycles=run_param.N_adapCycles;
    % post(ii).var_optim=run_param.var_optim;
    var_optim=adap_var_array(ii); % place holder since it is not saved
    post(ii).var_optim=var_optim;
    n_Dim=run_param.n_Dim;
    k_Order=run_param.k_Order;
    R_support=run_param.R_support;
    w_scale=run_param.w_scale;
    result_Var=run_param.result_Var;
    
    valid_ref=find(data_ref(:,1));
    points=data_ref(valid_ref,desVar); % Points where 
    
    %% Computing the monomial combinations symbolic expression
    % Calculates all the monomial combinations for an n dimensional space for
    % the recovery of a polynomial hypersurface of degree k.
    


    [mono_Sym, X_sym]=monomialSym(n_Dim,k_Order);

    mono_char=char(mono_Sym);
    mono_char=regexprep(mono_char,'_\d_','\(:,$0\)');
    mono_char=regexprep(mono_char,'_','');
    mono_char=regexprep(mono_char,'matrix([','');
    mono_char=regexprep(mono_char,']])',' ]');
    mono_char=regexprep(mono_char,'*','.*');
    mono_char=regexprep(mono_char,'\^','.\^');
    mono_char=regexprep(mono_char,'[1,','[ones(length(X(:,1)),1), ');

    fileID=fopen('MonomialFunc.m','w');
    lineFunc{1}='function F=MonomialFunc(X)';
    lineFunc{2}=['F=',mono_char,';'];

    for jj = 1:length(lineFunc)
        fprintf(fileID,'%s\n',lineFunc{jj});

    end

    fclose(fileID);

    for jj=1:length(X_sym)

        Lap_poly{jj}=char(diff(mono_Sym,X_sym(jj),2));
        Lap_poly{jj}=regexprep(Lap_poly{jj},'_\d_','\(:,$0\)');
        Lap_poly{jj}=regexprep(Lap_poly{jj},'_','');
        Lap_poly{jj}=regexprep(Lap_poly{jj},'matrix([','');
        Lap_poly{jj}=regexprep(Lap_poly{jj},'])','');
        Lap_poly{jj}=regexprep(Lap_poly{jj},'*','.*');
        Lap_poly{jj}=regexprep(Lap_poly{jj},'\^','.\^');
        Lap_poly{jj}=regexprep(Lap_poly{jj},'\d,||\d]','ones(length(X(:,1)),1)*$0');
        Lap_poly{jj}=['LapF(:,',num2str(jj),')=',Lap_poly{jj},'*Gamma;'];
    end


    fileID=fopen('MonomialLaplacian.m','w');

    lineFunc2{1}='function LaplaF=MonomialLaplacian(X,Gamma)';

    for jj=1:length(Lap_poly)
        lineFunc2{jj+1}=Lap_poly{jj};
    end
    lineFunc2{jj+2}='LaplaF=sum(LapF,2);';


    for jj = 1:length(lineFunc2)
        fprintf(fileID,'%s\n',lineFunc2{jj});

    end

    fclose(fileID);
    
    %% Load LHS Run
    
    
    load([run_dir,'OLHSrun.mat'])
    
    surrogateDat=SurrogatePost(points,data(:,[1 desVar result_Var(var_optim)]),R_support,w_scale);
    post(ii).surrogateDat.LHS=surrogateDat;
    
    [m,n]=size(surrogateDat);
    RMSE=sqrt((surrogateDat(:,end)-data_ref(valid_ref,result_Var(var_optim))).^2);
    post(ii).RMSE.LHS=sum(RMSE)/m;
    RMSErun=[0,N_LHS+n_Dim^2;post(ii).RMSE.LHS,post(ii).RMSE.LHS];
    
    clear data
    
    
    %% Perform Computations for Adaptive runs
    
    load([run_dir,'FullRun.mat'])
    
    adaptiveRuns=find(data(:,1)==2);
    if ~isempty(adaptiveRuns)
    for kk=adaptiveRuns(1):adaptiveRuns(end)
        
        surrogateDat=SurrogatePost(points,data(1:kk,[1 desVar result_Var(var_optim)]),R_support,w_scale);
        post(ii).surrogateDat.run(kk).dat=surrogateDat;

        [m,n]=size(surrogateDat);
        SquaredE=(surrogateDat(:,end)-data_ref(valid_ref,result_Var(var_optim))).^2;
        RMSErun=[RMSErun,[kk;sqrt(sum(SquaredE)/m)]]; 
        
        
        
    end
    end
    post(ii).RMSE.FullRun=RMSErun;
    
 
    
    
    disp('Done')
    
end

   %% Plot
    
%     line_format={'k-','b-.','r-.','g-.','b:','r:','g:','b--','r--','g--'};
%    
%     
%     for kk=1:3
%         cases=find(adap_var_array==kk);
%         cases(cases==kk)=[];
%         cases=[kk,cases];
%         figure
%         for ii=1:length(cases)
%             semilogy(post(cases(ii)).RMSE.FullRun(1,:),post(cases(ii)).RMSE.FullRun(2,:),line_format{ii},'linewidth',2);
%             hold on;
%         end
%     box=axis;
%     box(1:2)=[0,200];
%     axis(box)
%     end
    
%% save 

runSAVE(post,root_dir,'ResultsPost','result_data')


























