function [multiplier, veclin]=ScaleLineVecToLinear(vec, ratioEnd, maxVal)
    if ~exist('maxVal','var'); maxVal=1; end;
    flag = true;
    kk=0;
    while flag
        targetDistrib = reshape(kk+numel(vec):-1:kk+1,size(vec));
        if numel(vec) > 1
            targetDistrib(end) = targetDistrib(end-1)*ratioEnd;
        end

        veclin = targetDistrib/sum(targetDistrib)*sum(vec);
        kk=kk+1;
        flag = veclin(1)<maxVal;
        ratioEnd = 1-(1-ratioEnd)*0.8;
        if kk>20
            flag = false;
            warning('hit max iterate of Linear vec scaling')
        end
    end
    if vec(1)<vec(end)
        veclin = flip(veclin);
    end
    multiplier = veclin./vec;
end