

xDat=rand([10,1]);yDat=rand([10,1]);
N=11;
x=linspace(min(xDat),max(xDat),500)';

y=((x*ones(1,N)).^(ones(size(x))*(1:N)))*(((xDat*ones(1,N)).^(ones(size(xDat))*(1:N)))\yDat);

figure,plot(xDat,yDat,'*b',x,y,'r-');