function [newFill, nnum]=testnewreseed(fill)
    
    % identify elements = 1
    fillElmReplace = abs(fill-1)<1e-10;
    nnum =sum(fillElmReplace);
    sumVal = sum(fill);
    % Randomly reseed the elements =1 with a number in the range [0.5-1]
    popreplace = rand([1,nnum])*0.5+0.5;
    fill(fillElmReplace)=popreplace;
    % identify active fill variables(fill above 0)
    newSumVal=sum(fill);
    sup0=fill>0;
    % Calculate a fill correction to maintain (sum(newFill)==sum(oldFill))
    perVolDelta=(sumVal-newSumVal)/sum(sup0);
    rebalDelta=fill(sup0)+(1-fill(sup0))...
        /(sum((1-fill(sup0))))*(sum(sup0))*perVolDelta;

    % Assign newFill
    newFill = zeros(size(fill));
    newFill(sup0)=rebalDelta;
    if(abs(sum(newFill)-sumVal)>eps)
        abs(sum(newFill)-sum(fill))
        error('Fills should match after this step')
    end
end