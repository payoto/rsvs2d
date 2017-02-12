curvFunc=@(pi,pip1,pim1,s1,s2)(-pi*(s1+s2)+pip1*s2+pim1*s1)/(s1^2*s2+s2^2*s1);

testPos=[1 0.75; 0.75 1 ; -0.75 1; -1 0.75; -0.75 -1; 0.75 -1 ; 1 0.75; 0.75 1];
ratioA=[0.1 0.3 0.5 0.8 1 2 4 10 20 50];
for jj=1:numel(ratioA)
    testPos2=testPos*ratioA(jj);
    for ii=2:size(testPos,1)-1
        pi1=testPos2(ii);
        pip1=testPos2(ii+1);
        pim1=testPos2(ii-1);
        curv(ii,1:2)=curvFunc(pi1,pip1,pim1,norm(pi1-pip1),norm(pi1-pim1));

    end
    aTotCurv(jj)=sum(sqrt(sum((curv.^2),2)));
end



plot(ratioA,sqrt(aTotCurv).*ratioA)