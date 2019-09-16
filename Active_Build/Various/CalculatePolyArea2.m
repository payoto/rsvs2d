function [A,A2,matCalc]=CalculatePolyArea2(points)
    
    pointsVec=points';
    pointsVec=pointsVec(:);
    %plot(points(:,1),points(:,2));
    n=length(points(:,1));
    centreMat=eye(2*n);
    centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]));
    
    [rotDif]=[0 -1 0 1; 1 0 -1 0];
    normMat=zeros(2*n);
    for ii=1:n-1
        normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
    end
    ii=n;
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
    normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
    A=0.25*(normMat*pointsVec)'*(centreMat*pointsVec);
    matCalc = normMat'*centreMat;
    A2=0.25*pointsVec'*matCalc*pointsVec;
end