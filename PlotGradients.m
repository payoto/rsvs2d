function []=PlotGradients(paramoptim)
    
    for ii=1:numel(paramoptim)
        supportOptim=paramoptim(ii).optim.supportOptim;
       
        figure('Name',['Gradients_',ExtractVariables({'optimCase'},paramoptim(ii))]...
            ,'Position',[ 100 150 1400 700])
        subplot(2,2,1)
         grads=vertcat(supportOptim.hist(:).gradfk);
        surf(log10(abs(grads)))
        ylabel('iteration')
        xlabel('design variable')
        title('gradients')
        grad
        view(0,90)
        subplot(2,2,2)
         grads=vertcat(supportOptim.hist(:).prevDir);
        surf(log10(abs(grads)))
        ylabel('iteration')
        xlabel('design variable')
        title('previous direction')
        view(0,90)
        subplot(2,2,3)
         grads=vertcat(supportOptim.hist(:).prevStep);
        surf(((grads)),(abs((grads))))
        ylabel('iteration')
        xlabel('design variable')
        title('previous direction')
        view(0,90)
        subplot(2,2,4)
        plot()
        
    end
    
    
    
end