function []=ThresholdTriathlon(t1,rep,imCorr,offX,offY,mult)
    [mFin]=HenerateMask(t1,rep);
    
    figure
    subplot(1,2,1)
    imshow(imCorr)
    imCorr=double(imCorr)/255;
    mFin(mFin~=0)=nan;
    imCorr(offY+1:offY+size(mFin,1),offX+1:offX+size(mFin,2),:)=(repmat(mFin,[1 1 3])+1)...
        .*imCorr(offY+1:offY+size(mFin,1),offX+1:offX+size(mFin,2),:);
    imCorr=inpaint_nans(imCorr);
    subplot(1,2,2)
    imshow(imCorr)
end



function [mFin]=HenerateMask(t1,n)
    
    cRight=imread('\\ads.bris.ac.uk\filestore\myfiles\studentUG13\ap1949\Desktop\triathlon\logo.png');
cRight=double(cRight);
cRight=double(sum(cRight>t1,3)>0);
newMask=cRight;
oldMask=cRight;
for kk=1:n
    
    newMask(1:end-1,:)=newMask(1:end-1,:)+oldMask(2:end,:);
    newMask(2:end,:)=oldMask(1:end-1,:)+newMask(2:end,:);
    newMask(:,1:end-1)=newMask(:,1:end-1)+oldMask(:,2:end);
    newMask(:,2:end)=oldMask(:,1:end-1)+newMask(:,2:end);
    oldMask=newMask;
end

mFin=ones(size(cRight))*0.5;
mFin(logical(newMask))=0;
mFin(logical(cRight))=1;
mFin=mFin(20:end-15,103:end-103);


figure
imshow(mFin)
mFin=-(mFin-0.5);
    
end


function []=ThresholdTriathlon2(t1,t2,offX,offY,threshold)
    
    
    t1Off=double(t1(1:end-offY,1:end-offX,:))/255;
    t2Off=double(t2(1+offY:end,1+offX:end,:))/255;
    
    
    
    
    tfin=(((t1Off+t2Off)/2.*(1-abs((t1Off-t2Off)))))>threshold;
    tfin2=(((t1Off+t2Off)/2-threshold*abs((t1Off-t2Off))))*threshold;
    
    imshow(double(tfin2));
    
end