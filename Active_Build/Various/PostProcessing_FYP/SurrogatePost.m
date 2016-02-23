function surrogateDat=SurrogatePost(points,data,R_support,w_scale)
% Produces the surrogate model designated by RBF and Poly for
% the N^Dim values specified in vector ranges
%   POINTS are the points at which the surrogate model must be evaluated
%   DATA is the reduced data set (flagVar desVar adapVar) that is to be
%        analysed.
%


%% RBF model generation
%% Generation of correllation and polynomial matrices
% Generation of the RBF Surrogate model, This step involves the generation
% the R (correlation matrix) and the F matrix of monomials

col_index=2:length(data(1,:))-1;

for ii=1:length(col_index)
    
    bmax=max(data(:,col_index(ii)));
    ranges(ii,2)=bmax;
    bmin=min(data(:,col_index(ii)));
    ranges(ii,1)=bmin;
    if bmax~=bmin
    data(:,col_index(ii))=(data(:,col_index(ii))-bmin)/(bmax-bmin);
    else
        data(:,col_index(ii))=1;
    end
end

[m,n]=size(points);
pointsUH=(points-(ones(m,1)*ranges(:,1)'))./(ones(m,1)*(ranges(:,2)-ranges(:,1))');


desVar=2:length(data(1,:))-1;
valid_data=find(data(:,1));

R=PhiRBF(data(valid_data,desVar),data(valid_data,desVar),R_support,w_scale);

F=MonomialFunc(data(valid_data,desVar));



%% Direct Method
%Inline functions for Gamma and Beta coefficient finding

% Gamma Inline function
Gamma_RBF=inline('inv(F''*inv(R)*F)*F''*inv(R)*y','F','R','y');
% Beta Inline function
Beta_RBF=inline('inv(R)*(y-F*Gamma)','R','y','F','Gamma');


    
    Gamma_init=Gamma_RBF(F,R,data(valid_data,end));
    Beta_init=Beta_RBF(R,data(valid_data,end),F,Gamma_init);




       
%% RBF and polynomial terms generations

R=PhiRBF(pointsUH,data(valid_data,2:end-1), R_support,w_scale);

F=MonomialFunc(pointsUH);
surrogateDat=[points,R*Beta_init+F*Gamma_init];%




