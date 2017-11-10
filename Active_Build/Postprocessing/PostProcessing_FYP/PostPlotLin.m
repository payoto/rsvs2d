close all
clc

%% Declaration 
result_dir='C:\Users\Alexandre\Documents\Travail\Uni\IXP\Simulations\XFOIL6.99\RunSets\2D_50pts\ResultsPost\';

load([result_dir,'result_data.mat'])


[m,n]=size(data);

if m*n==max(m,n)

    for ii=1:n

        adap_var_array(ii)=data(ii).var_optim;


    end

    %% Plot



        line_format={'k-','b-.','r-.','g-.','b:','r:','g:','b--','r--','g--'};


        for kk=1:3
            cases=find(adap_var_array==kk);
            cases(cases==kk)=[];
            cases=[kk,cases];
            figure
            for ii=1:length(cases)
                semilogy(data(cases(ii)).RMSE.FullRun(1,:),data(cases(ii)).RMSE.FullRun(2,:),line_format{ii},'linewidth',2);
                hold on;
            end
        box=axis;
        box(1:2)=[min(post(cases(1)).RMSE.FullRun(1,:)),max(post(cases(1)).RMSE.FullRun(1,:))];
        axis(box)
        end
    
else
    line_format={'k-','b-.','r-.','g-.','b:','r:','g:','b--','r--','g--'};
    
    for jj=1:m
        
        for ii=1:n

            var_optim(jj,ii)=data(jj,ii).var_optim;
            
            [~,col(jj,ii)]=size(data(jj,ii).RMSE.FullRun);
            
            
        end
       
    end
    
    check_size=max(col,[],2);
   
    for jj=1:m
        cellRMSE{jj}=0;
        kk=0;
        for ii=1:n
            if col(jj,ii)==check_size(jj)
            cellRMSE{jj}=cellRMSE{jj}+data(jj,ii).RMSE.FullRun;
            kk=kk+1;
            end
            
        end
        cellRMSE{jj}=cellRMSE{jj}/kk;
    end
    
        for kk=1:3
        cases=find(var_optim(1:end,1)==kk);
%         cases(cases==kk)=[];
%         cases=[kk,cases'];
        figure
        ymin=[];
        ymax=[];
        for ii=1:length(cases)
                plot(cellRMSE{cases(ii)}(1,:),cellRMSE{cases(ii)}(2,:),line_format{ii},'linewidth',2);
                hold on;
                if isempty(ymin)
                    ymin=min(cellRMSE{cases(ii)}(2,:));
                    ymax=cellRMSE{cases(ii)}(2,1);
                else
                    ymin=min(ymin,min(cellRMSE{cases(ii)}(2,:)));
                    ymax=max(ymax,cellRMSE{cases(ii)}(2,1));
                end
        end
        ymin=ymin*10^(-0.1);
        ymax=ymax*10^(0.2);
        box=axis;
        box(3:4)=[ymin,ymax];
        axis(box)
        
        end
    
    
end

%% save 

%runSAVE(post,root_dir,'ResultsPost','result_data',1)