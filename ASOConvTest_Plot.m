
function []=ASOConvTest_Plot(ASOConvstruct)
    % Fill empties will Nan
    
    for ii=1:numel(ASOConvstruct)
        fields=fieldnames(ASOConvstruct{ii});
        for jj=[1,3:8]
            for kk=1:numel(ASOConvstruct{ii})
                if isempty(ASOConvstruct{ii}(kk).(fields{jj}))
                    ASOConvstruct{ii}(kk).(fields{jj})=nan;
                end
            end
        end
        
    end
    
    % Plot residuals
    hi=1;
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
    
    hi=hi+1;
    h(hi)=figure('Name','Iteration Plots');
    ax(3)=subplot(2,1,1);
    hold on
    ax(4)=subplot(2,1,2);
    hold on
    for ii=1:numel(ASOConvstruct)
        plot(ax(3),[ASOConvstruct{ii}.directNIter],'o-')
        plot(ax(4),[ASOConvstruct{ii}.adjointNIter],'o-')
        
    end
    
    
    hi=hi+1;
    h(hi)=figure('Name','ASOrunplots');
    numASOrun=zeros(size(ASOConvstruct));
    for ii=1:numel(ASOConvstruct)
        
        numASOrun(ii)=sum([ASOConvstruct{ii}.isAsoRun],'omitnan');
    end
    plot(numASOrun,'o-')
    
    % Plot iteration histories and adjoint histories
    
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
    h(hi)=figure('Name','Res_Rho_ Iteration histories');
    ax(lll+5)=subplot(1,2,1);
    hold on
    ax(lll+6)=subplot(1,2,2);
    hold on
    
    hi=hi+1;
    h(hi)=figure('Name','Res_RhoE_ Iteration histories');
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