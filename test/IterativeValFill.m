function [fill,isConstr]=IterativeMeanFill(fill,desVarRange,constrVal,volVec,totvol)
    %{
    isConstr=true;
    ratio=constrVal/(sum(fillStart.*volVec)/totvol);
    kk=0;
    n=length(fill);
    while ratio~=constrVal && kk<=n+1;
        maxFill=max(desVarRange);
        
        fillBound=((fill*ratio)>=maxFill);
        fill(fillBound)=maxFill;
        fill(~fillBound)=fill(~fillBound)*ratio;
        ratio=constrVal/(sum(fill.*volVec)/totvol);
        kk=kk+1;
    end
    if kk>n+1
        isConstr=false;
    end
    %}
    isConstr=true;
    ratio=2;
    kk=0;
    n=length(fill);
    maxFill=max(desVarRange);
    while ratio>1 && kk<=n+1;
        fillBound=(fill>=maxFill);
        ratio=(constrVal*totvol-sum(volVec(fillBound).*fill(fillBound)))/sum(volVec(~fillBound).*fill(~fillBound));
        fill(~fillBound)=min(fill(~fillBound)*ratio,maxFill);
        kk=kk+1;
    end
    if kk>n+1
        isConstr=false;
    end
end