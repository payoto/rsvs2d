
function [loop]=ConstantArea_Busemann(xMin,xMax,A,M)
    
    h=A/2/(xMax-xMin);
    Dx=xMax-xMin;
    tanDel=2*h/Dx;
    points=[xMin,0;xMax,0;xMin+(xMax-xMin)/2,h];
    
    [np]=GoldenSection_maxTanDel(0,pi/2,1e-6,M);
    [B]=GoldenSection_FindB(0,np,1e-6,M,tanDel);
    
    H=tan(B)*Dx/2;
    offset=H+h;
    points2=points;
    points2(:,2)=-points2(:,2)+offset/2;
    
    points(:,2)=points(:,2)-offset/2;
    points2=points2(end:-1:1,:);
    loop.subdivision=[points;points2];
    loop.isccw=true;
    
end

function [tanDel]=CalcTanDel(M,B)
    
    tanDel=2/tan(B)*(M^2*(sin(B)^2)-1)/(M^2*(1.4+cos(2*B))+2);
    
end

function [np]=GoldenSection_maxTanDel(lb,hb,tol,M)
    
    gr=(sqrt(5)-1)/2;
    
    c=hb-gr * (hb-lb);
    d=lb+gr * (hb-lb);
    
    while abs(c-d)>tol
        
        fc=-CalcTanDel(M,c);
        fd=-CalcTanDel(M,d);
        
        if fc<fd
            hb=d;
            
        else
            lb=c;
        end
        
        c=hb-gr * (hb-lb);
        d=lb+gr * (hb-lb);
    end
    
    np=(hb+lb)/2;
    
end

function [np]=GoldenSection_FindB(lb,hb,tol,M,targ)
    
    gr=(sqrt(5)-1)/2;
    
    c=hb-gr * (hb-lb);
    d=lb+gr * (hb-lb);
    
    while abs(c-d)>tol
        
        fc=abs(targ-CalcTanDel(M,c));
        fd=abs(targ-CalcTanDel(M,d));
        
        if fc<fd
            hb=d;
            
        else
            lb=c;
        end
        
        c=hb-gr * (hb-lb);
        d=lb+gr * (hb-lb);
    end
    
    np=(hb+lb)/2;
    
end