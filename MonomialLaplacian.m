function LaplaF=MonomialLaplacian(X,Gamma)
LapF(:,1)=[ones(length(X(:,1)),1)*0, ones(length(X(:,1)),1)*0, ones(length(X(:,1)),1)*0]*Gamma;
LapF(:,2)=[ones(length(X(:,1)),1)*0, ones(length(X(:,1)),1)*0, ones(length(X(:,1)),1)*0]*Gamma;
LaplaF=sum(LapF,2);
