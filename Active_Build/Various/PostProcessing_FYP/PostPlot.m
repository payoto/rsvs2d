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
        title_cell={'Graph 1: Evolution of RMSE with number of samples for the lift Coefficient';
            'Graph 2: Evolution of RMSE with number of samples for the lift Coefficient';
            'Graph 3: Evolution of RMSE with number of samples for the Moment Coefficient'};

        for kk=1:3
            cases=find(adap_var_array==kk);
            cases(cases==kk)=[];
            cases=[kk,cases];
            figure
            for ii=1:length(cases)
                line(ii,kk)=semilogy(data(cases(ii)).RMSE.FullRun(1,:),...
                    data(cases(ii)).RMSE.FullRun(2,:),line_format{ii},'linewidth',2);
                hold on;
                
                end
            box=axis;
            box(1:2)=[min(post(cases(1)).RMSE.FullRun(1,:)),max(post(cases(1)).RMSE.FullRun(1,:))];
            axis(box)
            
        end
    
else
    line_format={'k-','b-.','r-.','g-.','b:','r:','g:','b--','r--','g--'};
    title_cell={'Graph 1: Evolution of RMSE with number of samples for the lift Coefficient';
            'Graph 2: Evolution of RMSE with number of samples for the lift Coefficient';
            'Graph 3: Evolution of RMSE with number of samples for the Moment Coefficient'};
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
    for jj=1:m
        cellstd{jj}=[];
        kk=1;
        for ii=1:n
            if col(jj,ii)==check_size(jj)
            cellstd{jj}=[cellstd{jj};data(jj,ii).RMSE.FullRun(2,:)];
            kk=kk+1;
            end
            
        end
        cellstd{jj}=[data(jj,1).RMSE.FullRun(1,:);std(cellstd{jj})];
    end
        for kk=1:3
        cases=find(var_optim(1:end,1)==kk);
%         cases(cases==kk)=[];
%         cases=[kk,cases'];
        figure
        ymin=[];
        ymax=[];
        for ii=1:length(cases)
                line(ii,kk)=semilogy(cellRMSE{cases(ii)}(1,:),cellRMSE{cases(ii)}(2,:)...
                    ,line_format{ii},'linewidth',2);
                cell_legend{ii,kk}=['LHS: ',num2str(data(cases(ii),1).N_LHS),'; Adapting: ',num2str(data(cases(ii),1).N_adapCycles)];
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
        legend(cell_legend(:,kk))
        ylabel('RMSE')
        xlabel('Number of samples')
        title(title_cell(kk))
        end
    for kk=1:3
        cases=find(var_optim(1:end,1)==kk);
%         cases(cases==kk)=[];
%         cases=[kk,cases'];
        figure
        ymin=[];
        ymax=[];
        for ii=1:length(cases)
                plot(ii,kk)=semilogy(cellstd{cases(ii)}(1,:),cellstd{cases(ii)}(2,:)...
                    ,line_format{ii},'linewidth',2);
                cell_legend{ii,kk}=['LHS: ',num2str(data(cases(ii),1).N_LHS),'; Adapting: ',num2str(data(cases(ii),1).N_adapCycles)];
                hold on;
%                 if isempty(ymin)
%                     ymin=min(cellRMSE{cases(ii)}(2,:));
%                     ymax=cellRMSE{cases(ii)}(2,1);
%                 else
%                     ymin=min(ymin,min(cellRMSE{cases(ii)}(2,:)));
%                     ymax=max(ymax,cellRMSE{cases(ii)}(2,1));
%                 end
        end
%         ymin=ymin*10^(-0.1);
%         ymax=ymax*10^(0.2);
%         box=axis;
%         box(3:4)=[ymin,ymax];
%         axis(box)
        legend(cell_legend(:,kk))
        ylabel('\sigma')
        xlabel('Number of samples')
        title(title_cell(kk))
        end
    
end

%% save 

%runSAVE(data,result_dir,'fig','result_data',1)