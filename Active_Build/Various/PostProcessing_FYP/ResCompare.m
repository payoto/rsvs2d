
%%
load bigdata
Re_valid=find(data(:,3)==3000000);
size(Re_valid)
figure(1)
hold on
plot3(data(Re_valid,4),data(Re_valid,5),data(Re_valid,6),'g*','linewidth',2)
figure(2)
hold on
plot3(data(Re_valid,4),data(Re_valid,5),data(Re_valid,7),'g*','linewidth',2)
figure(3)
hold on
plot3(data(Re_valid,4),data(Re_valid,5),data(Re_valid,8),'g*','linewidth',2)
edit ResCompare

%%

valid_data=find(dat(:,1));
R=PhiRBF(dat(valid_data,desVar),dat(valid_data,desVar),R_support,w_scale);
F=PolyRBFchar(dat(valid_data,desVar),mono_char);
for ii=1:length(result_Var)
   
    Gamma_init(:,ii)=LSMPoly(F,dat(valid_data,result_Var(ii)));
    Delta=dat(valid_data,result_Var(ii))-F*Gamma_init(:,ii);
    Beta_init(:,ii)=(R)\Delta;
end


for ii=1:n_Dim
    ranges(ii,1)=min(dat(valid_data,desVar(ii)));
    ranges(ii,2)=max(dat(valid_data,desVar(ii)));
end

[surDat,surStruct]=Surrogate(ranges,50,dat(valid_data,desVar),Beta_init,Gamma_init,mono_char);
disp('Done.')

t=[t,toc];
ResPlot(surStruct,X_plot,X_default,dat(valid_data,[desVar(X_plot),result_Var]))
