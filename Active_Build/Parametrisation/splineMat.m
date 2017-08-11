%% Recursive basis function generator
function Z = SplineBasis(i,p,u,U)  % from dom
    i = i+1;
    if p==0
        % Z=zeros(1,length(u));
        for j=1:length(u);
            if j==length(u)
                if u(j)>=U(i) && u(j)<=U(i+1)
                    Z(j)=1;
                else
                    Z(j) = 0;
                end
            else
                if u(j)>=U(i) && u(j)<U(i+1)
                    Z(j)=1;
                else
                    Z(j)=0;
                end
            end
        end
    else
        %         Z1=((u-U(i))/(U(i+p)-U(i))).*N(i-1,p-1,u,U);
        %         Z2=((U(i+p+1)-u)/(U(i+p+1)-U(i+1))).*N((i+1)-1,p-1,u,U);
        %         if isnan(Z1)==1;
        %             if isnan(Z2)==1;
        %                 Z=zeros(size(u));
        %             else
        %                 Z=Z2;
        %             end
        %         else
        %             if isnan(Z2)==1;
        %                 Z=Z1;
        %             else
        %                 Z=Z1+Z2;
        %             end
        %         end
        
        t1n = (u-U(i));
        t1d = (U(i+p)-U(i));
        t2n = (U(i+p+1)-u);
        t2d = (U(i+p+1)-U(i+1));
        
        b1 = N(i-1,p-1,u,U);
        b2 = N((i+1)-1,p-1,u,U);
        
        if t1d ~= 0
            t1 = t1n/t1d;
        else
            t1 = 0;
        end
        
        if t2d ~= 0
            t2 = t2n/t2d;
        else
            t2 = 0;
        end
        
        Z = t1 .* b1 + t2 .* b2;
        
        
        
    end
end


function Nmat = splineMat(dataKnots, splineKnots, order)
    
    splineKnots = [zeros(1,order) , splineKnots, ones(1,order)];
    L = length(splineKnots) - (1+ order);
    
    Nmat = zeros(length(dataKnots),L);
    for i=1:L
        Nmat(:,i) = N(i-1,order,dataKnots, splineKnots);
    end
    
end
