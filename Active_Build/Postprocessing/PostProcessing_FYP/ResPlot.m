function []=ResPlot(surStruct,dims,static,data)
% Produces a 3D plot for the surrogate model designated by RBF and Poly for
% the N^Dim values specified in vector ranges
%   RANGES: is a D*2 matrix containing the lower and upper bounds of each
%           variable
%        N: is the number of graduations in each dimension
% MODELDAT: is the matrix containign the data points used to compute the
%           RBF model
%      RBF: is the coefficients associated with the RBF model
%  polySym: Is the symbolic vector of monomials
%     poly: is the set of coefficients for the polynomial



for ii=1:length(static)
    compare=abs(surStruct.X(ii).dat-ones(size(surStruct.X(ii).dat))*static(ii));
    [~,I_static(ii)]=min(compare);
end

if length(dims)>1
    dims=sort(dims);
    [~,n]=size(surStruct.Grid);

    [X,Y]=meshgrid(surStruct.X(dims(1)).dat,surStruct.X(dims(2)).dat);

    m=size(surStruct.Grid(1).res);
    k=length(m);
    m=m(1);


    for ii=1:k
        if sum(dims==ii)
            Index{ii}=1:m;
        else
            Index{ii}=I_static(ii);
        end
    end

    orig=find(data(:,end)==1);
    adapt=find(data(:,end)==2);

    for ii=1:n
        Z=surStruct.Grid(ii).res(Index{:});
        Z=squeeze(Z);
        figure
        surf(X,Y,Z')
        hold on

        if length(static)==2
            plot3(data(orig,1),data(orig,2),data(orig,2+ii),'k*','linewidth',3)
            hold on
            plot3(data(adapt,1),data(adapt,2),data(adapt,2+ii),'w*','linewidth',3)
            hold on
        end
    end
else
    figure
   plot(surStruct.X.dat,surStruct.Grid.res(:,1))
   hold on
   plot(data(:,1),data(:,2),'k*')
    axis([0 1 -.2 1.2])
end



