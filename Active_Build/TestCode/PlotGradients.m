function [h]=PlotGradients(paramoptim,optimstruct)
    
    normVec=@(vec) sqrt(sum(vec.^2,2));
    normaliseArray=@(array)  array./(sqrt(sum(array.^2,2))*ones(1,size(array,2))); 
    
    supportOptim=paramoptim.optim.supportOptim;
    
    h=figure('Name',['Gradients_',ExtractVariables({'optimCase'},paramoptim)]...
        ,'Position',[ 100 150 1400 700]);
    symDesVarList=ExtractVariables({'symDesVarList'},paramoptim);
    rmCol=symDesVarList(2,:);
    subplot(2,2,1)
    grads=vertcat(supportOptim.hist(:).gradfk);
    grads(:,rmCol)=[];
    surf(log10(abs(grads)))
    ylabel('iteration')
    xlabel('design variable')
    title('gradients')
    
    view(0,90)
    subplot(2,2,2)
    grads=vertcat(supportOptim.hist(:).prevDir);
    grads(:,rmCol)=[];
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
    hold on
    grads=vertcat(supportOptim.hist(:).gradfk);
    grads(:,rmCol)=[];
    gradsm1=vertcat(supportOptim.hist(:).gradfkm1);
    gradsm1(:,rmCol)=[];
    directionChange.Grad=dot(normaliseArray(grads),normaliseArray(gradsm1),2);
    directionChange.Grad(1)=1;
    gradStep=vertcat(supportOptim.hist(:).prevStep);
    gradStep(:,rmCol)=[];
    directionChange.Step=dot(normaliseArray(grads),normaliseArray(gradStep),2);
    gradDir=vertcat(supportOptim.hist(:).prevDir);
    gradDir(:,rmCol)=[];
    directionChange.Dir=dot(normaliseArray(grads),normaliseArray(gradDir),2);
    
    kk=1;
    fillPrec=zeros(size(optimstruct(1).population(1).fill));
    for ii=1:2:numel(optimstruct)
        changePos(kk)=normVec(fillPrec-optimstruct(ii).population(1).fill);
        fillInf(kk,1:numel(optimstruct(ii).population(1).fill))...
            =optimstruct(ii).population(1).fill;
        fillPrec=optimstruct(ii).population(1).fill;
        kk=kk+1;
    end
    
    
    
    step=1:numel(supportOptim.hist);
    ii=1;
    l(ii)=plot(step,directionChange.Grad);
    l(ii).DisplayName='Change in direction - gradient';
    ii=ii+1;
    l(ii)=plot(step,directionChange.Step);
    l(ii).DisplayName='gradient // direction';
    ii=ii+1;
    l(ii)=plot(step,directionChange.Dir);
    l(ii).DisplayName='gradient // step';
    ii=ii+1;
    try
        scaleEvol=[supportOptim.hist(:).scale];
        l(ii)=plot(step,scaleEvol);
        l(ii).DisplayName='scale';
        ii=ii+1;
    catch
    end
    try
        iterEvol=[supportOptim.hist(:).iter];
        l(ii)=plot(step,iterEvol);
        l(ii).DisplayName='iter since last refresh';
        ii=ii+1;
    catch
    end
    l(ii)=plot(step,changePos);
    l(ii).DisplayName='Length of movement';
    ii=ii+1;
    legend(l)
    
end