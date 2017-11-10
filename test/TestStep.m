function [stepLength]=TestStep(f,worker,direction,test)
    
    
    if test==1
        
    
    
    
    [~,bestPoint]=min(f);
    stepLengths=1./2.^[inf,(worker-2):-1:0];
    
    if bestPoint==1
        bestPoint=2;
    end
    if bestPoint==worker
        bestPoint=worker-1;
    end
    prevPoint=bestPoint-1;
    nextPoint=bestPoint+1;
    
    indSearch=[prevPoint,bestPoint,nextPoint];
    
    
    %pp=spline(stepLengths(indSearch),f(indSearch));
    [pp]=ParabolicFit(stepLengths(indSearch),f(indSearch));
    iTest=linspace(stepLengths(indSearch(1)),stepLengths(indSearch(3)),1000);
    switch direction
        case 'min'
            [targObj,indexLoc]=min(ParabolicVal(pp,iTest));
        case 'max'
            [targObj,indexLoc]=max(ParabolicVal(pp,iTest));
    end
    stepLength=iTest(indexLoc);
    
    if stepLength==0
        warning('Step Length is stagnant this iteration')
    end
    
    elseif test==2
        
        [~,bestPoint]=min(f);
        stepLengths=1./2.^[inf,(worker-2):-1:0];
        
        if bestPoint==1
            bestPoint=2;
        end
        if bestPoint==worker
            bestPoint=worker-1;
        end
        prevPoint=bestPoint-1;
        nextPoint=bestPoint+1;
        
        indSearch=[prevPoint,nextPoint,prevPoint];
        
        
        pp=spline(stepLengths,f);
        
        iTest=linspace(0,1,1000);
        switch direction
            case 'min'
                [targObj,indexLoc]=min(ppval(pp,iTest));
            case 'max'
                [targObj,indexLoc]=max(ppval(pp,iTest));
        end
        stepLength=iTest(indexLoc);
        
        if stepLength==0
            warning('Step Length is stagnant this iteration')
        end
        
    end
    
    
end

function [coeff]=ParabolicFit(xI,yI)
    
    parabola=@(x) [x.^2, x ,ones(size(x))];
    
    xI=reshape(xI,[numel(xI),1]);
    yI=reshape(yI,[numel(yI),1]);
    R=parabola(xI);
    
    if numel(xI)==3
        coeff=R\yI;
    else
        warning('Least Square fit used')
        coeff=(R'*R)\R'*yI;
        
    end
    
end

function [yy]=ParabolicVal(coeff,xx)
    parabola=@(x) [x.^2, x ,ones(size(x))];
    
    xx=reshape(xx,[numel(xx),1]);
    
    yy=parabola(xx)*coeff;
end