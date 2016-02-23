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

root_dir='Z:\Documents\Work\year4\research\FYP\MATLAB\Results\XFOIL6.99\RunSets4\2D_50pts\';
data_file='bigdata2D2.mat';
% need to specify which variables have been otimised
% 
adap_var_array_LHS=[1, 2 ,3];



d = dir(root_dir);
isub = [d(:).isdir]; %# returns logical vector
nameSet = {d(isub).name}';
logiSet=regexp(nameSet,'set');



for ii=1:length(nameSet)
    if isempty(logiSet{ii})
        logiSet{ii}=0;
    end
end
logiSet=logical(cell2mat(logiSet));

nameSet=nameSet(logiSet);
% numFolds=regexprep(nameFolds,'run','');
% numFolds=cell2mat(numFolds);


n_set=sum(logiSet);


% Run folder to be analyzed
set_dir=[root_dir,'set1\'];
d = dir(set_dir);
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

load([root_dir,data_file])

data_ref=bigdata;

clear bigdata

%% Analysis of Subfolders



for mm=1:n_set
    
    % Run folder to be analyzed
    set_dir=[root_dir,'set',num2str(mm),'\'];
    d = dir(set_dir);
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
    
    col_num=1;
for ii=1:n_run
    if ii==1
        %to be repeated for all optim_var
        repeat=length(adap_var_array_LHS);
        
    else
        repeat=1;
        
    end
    
    for ii_repeat=1:repeat
        
    
    
        %% Subdirectory specific variable extraction

        fprintf('Run %i out of %i started\n',ii,n_run)

        run_dir=[root_dir,nameSet{mm},'\',nameFolds{ii},'\'];

        load([run_dir,'RunParam.mat'])

        run_param=data;
        clear data

        post(col_num,mm).desVar=run_param.desVar;
        desVar=run_param.desVar;
        post(col_num,mm).des_constant=run_param.des_constant;
        des_constant=run_param.des_constant;
        post(col_num,mm).N_LHS=run_param.N_LHS;
        N_LHS=run_param.N_LHS;
        post(col_num,mm).N_adapCycles=run_param.N_adapCycles;
        N_adapCycles=run_param.N_adapCycles;
        n_Dim=run_param.n_Dim;
        k_Order=run_param.k_Order;
        R_support=run_param.R_support;
        w_scale=run_param.w_scale;
        result_Var=run_param.result_Var;
        
        if ii==1
            var_optim=adap_var_array_LHS(ii_repeat);
            post(col_num,mm).var_optim=var_optim;
        else
            post(col_num,mm).var_optim=run_param.var_optim;
            %var_optim=adap_var_array(ii); % place holder since it is not saved
            var_optim=post(col_num,mm).var_optim;
        end
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
        post(col_num,mm).surrogateDat.LHS=surrogateDat;

        [m,n]=size(surrogateDat);
        RMSE=((surrogateDat(:,end)-data_ref(valid_ref,result_Var(var_optim))).^2);
        post(col_num,mm).RMSE.LHS=sqrt(sum(RMSE)/m);
        RMSErun=[0,N_LHS+n_Dim^2;post(col_num,mm).RMSE.LHS,post(col_num,mm).RMSE.LHS];

        clear data


        %% Perform Computations for Adaptive runs

        load([run_dir,'FullRun.mat'])

     max_ind=max(find(data(:,1)>0));
        for kk=(N_LHS+n_Dim^2)+1:max_ind

            surrogateDat=SurrogatePost(points,data(1:kk,[1 desVar result_Var(var_optim)]),R_support,w_scale);
            post(col_num,mm).surrogateDat.run(kk).dat=surrogateDat;

            [m,n]=size(surrogateDat);
            SquaredE=(surrogateDat(:,end)-data_ref(valid_ref,result_Var(var_optim))).^2;
            RMSErun=[RMSErun,[kk;sqrt(sum(SquaredE)/m)]]; 



        end

        post(col_num,mm).RMSE.FullRun=RMSErun;


    col_num=col_num+1;
    end
        
    
end
end

%% Check
[m,n]=size(post);
for jj=1:m
        
        for ii=1:n

            var_optim(jj,ii)=post(jj,ii).var_optim;
            
            [~,col(jj,ii)]=size(post(jj,ii).RMSE.FullRun);
            
            
        end
       
    end
col

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
%             semilogy(post(cases(ii,mm)).RMSE.FullRun(1,:),post(cases(ii,mm)).RMSE.FullRun(2,:),line_format{ii},'linewidth',2);
%             hold on;
%         end
%     box=axis;
%     box(1:2)=[0,200];
%     axis(box)
%     end
    
%% save 

runSAVE(post,root_dir,'ResultsPost','result_data')


























