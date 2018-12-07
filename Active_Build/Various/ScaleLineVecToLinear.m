function [multiplier, veclin]=ScaleLineVecToLinear(vec, ratioEnd)
    
   sum(vec);
   targetDistrib = (numel(vec):-1:1);
   targetDistrib(end) = ratioEnd;
    
   veclin = targetDistrib/sum(targetDistrib)*sum(vec);
   
end