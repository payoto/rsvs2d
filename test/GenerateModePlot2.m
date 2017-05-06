
for ii=2:3;
    h=findobj(ii,'type','axes');
    if numel(h)>0
        for jj=1:numel(h)
            h(jj).Position=hRef(3).Position;
            h(jj).OuterPosition=hRef(3).OuterPosition;
        end
    end
    %h=findobj(ii);
    %h(1).Position=hRef(1).Position;
end
