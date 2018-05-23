
%%
fig=5
l=findobj(fig,'type','line');
leg=findobj(fig,'type','legend');
ax=findobj(fig,'type','axes')
lineLeg=findobj(ax(3),'type','line')
legend(flip(lineLeg))
leg=findobj(fig,'type','legend');


%%
tabReplace={{'d','-','none','-'};{'+','-','none','--'};{'.','-','.','--'}};
for ii=4:numel(l)
    flag=true;
    for jj=1:numel(tabReplace)
        if flag
            if l(ii).LineStyle==tabReplace{jj}{2} && l(ii).Marker(1)==tabReplace{jj}{1}
                l(ii).LineStyle=tabReplace{jj}{4};
                l(ii).Marker=tabReplace{jj}{3};
                flag=false;
            end
        end
    end
end

dispName=regexprep({lineLeg.DisplayName},'^.*conv ([0-9]) 7 ([0-9]*)','M-L: $1 to 7 (iter lim: $2)');
dispName=regexprep(dispName,'^.*conv ([0-9]) 1 ([0-9]*)','Level: $1 (iter lim: 100)');
[lineLeg.DisplayName]=deal(dispName{:})

ax(2).Position(1)=0.572
ax(3).Position(1)=0.572
ax(3).Position(1)=0.1
ax(1).Position(1)=0.572
ax(1).Position(1)=0.60
ax(2).Position(1)=0.60
ax(2).Position(3)=0.38
ax(1).Position(3)=0.38
ax(3).YAxis.Label.String='$C_D$'
ax(3).XAxis.Label.String='$SNOPT Major Iteration$'
ax(2).Position(2)=ax(1).Position(2)
ax(2).Position(2)=ax(3).Position(2)

FigureTextTools([fig],'fontSize',12,'interpreter','latex')
QuickFigSave([fig])