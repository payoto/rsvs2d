function [out]=ChokeBoundModelFunction(x)
    
    out = 0.5*sum(x,2) + sum(x.^2,2) + (sum(x.^2,2)<1).*(10- 3*sum(x.^2,2));