
function []=ASOConvTest_Plot(ASOConvstruct,isIter)
    % Fill empties will Nan
    
    if nargin<2
        isIter=false;
    end
    for ii=1:numel(ASOConvstruct)
        fields=fieldnames(ASOConvstruct{ii});
        if(~isfield(ASOConvstruct{ii},'directResidual'))
            for kk=1:numel(ASOConvstruct{ii})
                ASOConvstruct{ii}(kk).directResidual=ASOConvstruct{ii}(kk).residual.direct;
            end
        end
        if(~isfield(ASOConvstruct{ii},'adjointResidual'))
            for kk=1:numel(ASOConvstruct{ii})
                ASOConvstruct{ii}(kk).adjointResidual=ASOConvstruct{ii}(kk).residual.adjoint;
            end
        end
        for jj=[1,3:8]
            for kk=1:numel(ASOConvstruct{ii})
                try
                    ASOConvstruct{ii}(kk).directMaxRes=ASOConvstruct{ii}(kk).directIterHist.maxResidual(end);
                catch
                    ASOConvstruct{ii}(kk).directMaxRes=[];
                end
                try
                    ASOConvstruct{ii}(kk).adjointMaxRes=ASOConvstruct{ii}(kk).adjointIterHist.maxResidual(end);
                catch
                    ASOConvstruct{ii}(kk).adjointMaxRes=[];
                end
                try
                    ASOConvstruct{ii}(kk).directT=sum(ASOConvstruct{ii}(kk).directIterHist.Time_s_);
                catch
                    ASOConvstruct{ii}(kk).directT=[];
                end
                try
                    ASOConvstruct{ii}(kk).adjointT=sum(ASOConvstruct{ii}(kk).adjointIterHist.Time_s_);
                catch
                    ASOConvstruct{ii}(kk).adjointT=[];
                end
                
                if isempty(ASOConvstruct{ii}(kk).(fields{jj}))
                    ASOConvstruct{ii}(kk).(fields{jj})=nan;
                end
            end
        end
        
    end
    
    % Plot residuals
    hi=1;
    try
        h(hi)=figure('Name','Residual Plots');
        ax(1)=subplot(1,2,1);
        cOrd=get(ax(1),'colororder');
        hold on
        ax(2)=subplot(1,2,2);
        hold on
        for ii=1:numel(ASOConvstruct)
            plot(ax(1),[ASOConvstruct{ii}.directResidual],'o-')
            plot(ax(2),[ASOConvstruct{ii}.adjointResidual],'o-')
            
        end
    catch
        warning('Residual plots failed')
    end
    
    hi=hi+1;
    try
        h(hi)=figure('Name','Iteration Plots');
        ax(3)=subplot(2,1,1);
        hold on
        ax(4)=subplot(2,1,2);
        hold on
        for ii=1:numel(ASOConvstruct)
            plot(ax(3),[ASOConvstruct{ii}.directNIter],'o-')
            plot(ax(4),[ASOConvstruct{ii}.adjointNIter],'o-')
            
        end
    catch
        warning('Iteration plots failed')
    end
    
    hi=hi+1;
    h(hi)=figure('Name','Max residual Plots');
    ax(3)=subplot(2,1,1);
    hold on
    ax(4)=subplot(2,1,2);
    hold on
    for ii=1:numel(ASOConvstruct)
        plot(ax(3),[ASOConvstruct{ii}.directMaxRes],'o-')
        plot(ax(4),[ASOConvstruct{ii}.adjointMaxRes],'o-')
        
    end
    hi=hi+1;
    try
        h(hi)=figure('Name','Computational Expense plots');
        ax(3)=subplot(2,1,1);
        hold on
        ax(4)=subplot(2,1,2);
        hold on
        for ii=1:numel(ASOConvstruct)
            l=plot(ax(3),[ASOConvstruct{ii}.directT],'o-');
            plot(ax(3),[1,numel(ASOConvstruct{ii})],mean([ASOConvstruct{ii}.directT])*[1 1],'--','color',l.Color)
            l=plot(ax(4),[ASOConvstruct{ii}.adjointT],'o-');
            plot(ax(4),[1,numel(ASOConvstruct{ii})],mean([ASOConvstruct{ii}.adjointT])*[1 1],'--','color',l.Color)
            compExpense(1,ii)=mean([ASOConvstruct{ii}.directT]);
            compExpense(2,ii)=mean([ASOConvstruct{ii}.adjointT]);
            compExpense(3,ii)=compExpense(1,ii)+compExpense(2,ii);
            
        end
    catch
        warning('Computational expense failed')
    end
    
    hi=hi+1;
    try
        h(hi)=figure('Name','ASOrunplots');
        subplot(2,1,1)
        numASOrun=zeros(size(ASOConvstruct));
        for ii=1:numel(ASOConvstruct)
            
            numASOrun(ii)=sum([ASOConvstruct{ii}.isAsoRun],'omitnan');
            
        end
        plot(numASOrun,'o-')
        ylabel('Number of ASO runs')
        subplot(2,1,2)
        plot(compExpense','o-')
        ylabel('Computational Time (s)')
        
    catch
        warning('ASORunplots plots failed')
    end
    
    % Plot iteration histories and adjoint histories
    
    if isIter
        hi=hi+1;
        h(hi)=figure('Name','Res_Rho_ Iteration histories');
        ax(5)=subplot(1,2,1);
        hold on
        ax(6)=subplot(1,2,2);
        hold on
        
        hi=hi+1;
        h(hi)=figure('Name','Res_RhoE_ Iteration histories');
        ax(7)=subplot(1,2,1);
        hold on
        ax(8)=subplot(1,2,2);
        hold on
        for ii=1:numel(ASOConvstruct)
            
            for jj=1:numel(ASOConvstruct{ii})
                try
                    plot3(ax(5),ASOConvstruct{ii}(jj).directIterHist.Iter,...
                        jj*ones(size(ASOConvstruct{ii}(jj).directIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).directIterHist.Res_Rho_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                try
                    plot3(ax(6),ASOConvstruct{ii}(jj).adjointIterHist.Iter,...
                        jj*ones(size(ASOConvstruct{ii}(jj).adjointIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).adjointIterHist.Res_Psi_Rho_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                try
                    plot3(ax(7),ASOConvstruct{ii}(jj).directIterHist.Iter,...
                        jj*ones(size(ASOConvstruct{ii}(jj).directIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).directIterHist.Res_RhoE_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                try
                    plot3(ax(8),ASOConvstruct{ii}(jj).adjointIterHist.Iter,...
                        jj*ones(size(ASOConvstruct{ii}(jj).adjointIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).adjointIterHist.Res_Psi_E_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                
            end
            
        end
        for ii=5:8
            box=axis(ax(ii));
            box(5:6)=[-10 5];
            box(1:2)=[0 3000];
            axis(ax(ii),box)
        end
        
        lll=4;
        hi=hi+1;
        h(hi)=figure('Name','Res_Rho_ Time histories');
        ax(lll+5)=subplot(1,2,1);
        hold on
        ax(lll+6)=subplot(1,2,2);
        hold on
        
        hi=hi+1;
        h(hi)=figure('Name','Res_RhoE_ Time histories');
        ax(lll+7)=subplot(1,2,1);
        hold on
        ax(lll+8)=subplot(1,2,2);
        hold on
        for ii=1:numel(ASOConvstruct)
            
            for jj=1:numel(ASOConvstruct{ii})
                try
                    plot3(ax(lll+5),cumsum(ASOConvstruct{ii}(jj).directIterHist.Time_s_),...
                        jj*ones(size(ASOConvstruct{ii}(jj).directIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).directIterHist.Res_Rho_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                try
                    plot3(ax(lll+6),cumsum(ASOConvstruct{ii}(jj).adjointIterHist.Time_s_),...
                        jj*ones(size(ASOConvstruct{ii}(jj).adjointIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).adjointIterHist.Res_Psi_Rho_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                try
                    plot3(ax(lll+7),cumsum(ASOConvstruct{ii}(jj).directIterHist.Time_s_),...
                        jj*ones(size(ASOConvstruct{ii}(jj).directIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).directIterHist.Res_RhoE_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                try
                    plot3(ax(lll+8),cumsum(ASOConvstruct{ii}(jj).adjointIterHist.Time_s_),...
                        jj*ones(size(ASOConvstruct{ii}(jj).adjointIterHist.Iter)),...
                        ASOConvstruct{ii}(jj).adjointIterHist.Res_Psi_E_,'-','color',...
                        cOrd(mod(ii-1,size(cOrd,1))+1,:));
                catch
                end
                
            end
            
        end
        for ii=5:8
            box=axis(ax(lll+ii));
            box(5:6)=[-10 5];
            box(1:2)=[0 100];
            axis(ax(lll+ii),box)
        end
    end
    
end