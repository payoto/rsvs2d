%% Code to generate a 3D plot of function

function []=Test_ProblemDerivatives(derivtenscalc,numStep,stepSize,eps)
    if ~exist('eps','var'), eps=0;end
    jj=0:numStep-1;
    dVec=1-jj*stepSize;
    [d_i,d_m]=meshgrid(dVec,dVec);
    for jj=1:numStep
        for kk=1:numStep
            derivtenscalc.d_i=dVec(jj)-eps;
            derivtenscalc.d_m=dVec(kk)-eps;
            derivtenscalc.normFi=sqrt(sum(((derivtenscalc.g1_i+derivtenscalc.Dg_i*derivtenscalc.d_i)-(derivtenscalc.g1_m+derivtenscalc.Dg_m*derivtenscalc.d_m)).^2));
            [derivtenscalc3(jj,kk)]=CalculateDerivatives(derivtenscalc);
        end
    end
    
    di=zeros(numStep);
    dm=zeros(numStep);
    Fi=zeros(numStep);
    dfiddi=zeros(numStep);
    dfiddm=zeros(numStep);
    d2fiddi2=zeros(numStep);
    d2fiddm2=zeros(numStep);
    d2fiddim=zeros(numStep);
    
    for jj=1:numStep
        di(jj,:)=[derivtenscalc3(jj,:).d_i];
        dm(jj,:)=[derivtenscalc3(jj,:).d_m];
        Fi(jj,:)=[derivtenscalc3(jj,:).normFi];
        dfiddi(jj,:)=[derivtenscalc3(jj,:).dfiddi];
        dfiddm(jj,:)=[derivtenscalc3(jj,:).dfiddm];
        
        d2fiddi2(jj,:)=[derivtenscalc3(jj,:).d2fiddi2];
        d2fiddm2(jj,:)=[derivtenscalc3(jj,:).d2fiddm2];
        d2fiddim(jj,:)=[derivtenscalc3(jj,:).d2fiddim];
    end
    for jj=2:numStep-1
        for kk=2:numStep-1
%             dfiddiFD(jj-1,kk-1)=(Fi(jj+1,kk)-Fi(jj-1,kk))/(di(jj+1,kk)-di(jj-1,kk));
%             dfiddmFD(jj-1,kk-1)=(Fi(jj,kk+1)-Fi(jj,kk-1))/(dm(jj,kk+1)-dm(jj,kk-1));
            d2fiddi2FD(jj-1,kk-1)=(dfiddi(jj+1,kk)-dfiddi(jj-1,kk))/(di(jj+1,kk)-di(jj-1,kk));
            d2fiddm2FD(jj-1,kk-1)=(dfiddm(jj,kk+1)-dfiddm(jj,kk-1))/(dm(jj,kk+1)-dm(jj,kk-1));
            d2fiddmiFD(jj-1,kk-1)=(dfiddi(jj,kk+1)-dfiddi(jj,kk-1))/(dm(jj,kk+1)-dm(jj,kk-1));
%             d2fiddi2FD(jj,kk)=
%             d2fiddm2FD(jj,kk)=
%             d2fiddimFD(jj,kk)=
        end
    end
    figure
    surf(1-di,1-dm,Fi,'linestyle','none')
    title('Norm (Fi)')
    xlabel('di')
    ylabel('dm')
    
    
    figure
    subplot(1,2,1)
    surf(1-di,1-dm,dfiddi,'linestyle','none')
    title('1st di derivative')
    xlabel('di')
    ylabel('dm')
    subplot(1,2,2)
    surf(1-di,1-dm,dfiddm,'linestyle','none')
    title('1st dm derivative')
    xlabel('di')
    ylabel('dm')
    
    
    figure
    subplot(1,3,1)
    [h(1)]=SurfFor2ndDeriv(1-di,1-dm,d2fiddi2);
    hold on
    %surf(1-di(2:end-1,2:end-1),1-dm(2:end-1,2:end-1),d2fiddi2FD,log10(abs(d2fiddi2FD)),'linestyle','none')
    title('2nd di derivative')
    xlabel('di')
    ylabel('dm')
    subplot(1,3,3)
    [h(2)]=SurfFor2ndDeriv(1-di,1-dm,d2fiddm2);
    hold on
    %surf(1-di(2:end-1,2:end-1),1-dm(2:end-1,2:end-1),d2fiddm2FD,log10(abs(d2fiddm2FD)),'linestyle','none')
    title('2nd dm derivative')
    xlabel('di')
    ylabel('dm')
    subplot(1,3,2)
    [h(3)]=SurfFor2ndDeriv(1-di,1-dm,d2fiddim);
    title('dmdi derivative')
    xlabel('di')
    ylabel('dm')
    hold on
    %surf(1-di(2:end-1,2:end-1),1-dm(2:end-1,2:end-1),d2fiddmiFD,log10(abs(d2fiddmiFD)),'linestyle','none')
    
end

function[]=vjkfd()
    for jj=1:1000
        conv=-1e-6-jj*10^-8;
        derivtenscalc(ii).d_i=1+conv;
        derivtenscalc(ii).d_m=1;
        derivtenscalc(ii).normFi=sqrt(sum(((derivtenscalc(ii).g1_i+derivtenscalc(ii).Dg_i*derivtenscalc(ii).d_i)-(derivtenscalc(ii).g1_m+derivtenscalc(ii).Dg_m*derivtenscalc(ii).d_m)).^2));
        [derivtenscalc(ii).a_i,...
            derivtenscalc(ii).a_m,...
            derivtenscalc(ii).a_im,...
            derivtenscalc(ii).b_i,...
            derivtenscalc(ii).b_m,...
            derivtenscalc(ii).c]=...
            Calc_LengthDerivCoeff(...
            derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m,...
            derivtenscalc(ii).g1_i,derivtenscalc(ii).g1_m);
        [derivtenscalc3(jj)]=CalculateDerivatives(derivtenscalc(ii));
    end
    di=[derivtenscalc3(:).d_i];
    dm=[derivtenscalc3(:).d_m];
    Fi=[derivtenscalc3(:).normFi];
    dfiddi=[derivtenscalc3(:).dfiddi];
    dfiddm=[derivtenscalc3(:).dfiddm];
    
    d2fiddi2=[derivtenscalc3(:).d2fiddi2];
    d2fiddm2=[derivtenscalc3(:).d2fiddm2];
    d2fiddim=[derivtenscalc3(:).d2fiddim];
    for jj=1:length(di)
        ip=min([length(di),jj+1]);
        im=max([1,jj-1]);
        dfiddiFD(jj)=(Fi(ip)-Fi(im))/(di(ip)-di(im));
    end
    for jj=2:length(di)-1
        dfiddi2FD(jj-1)=(Fi(jj-1)+Fi(jj+1)-2*Fi(jj))/(di(jj)-di(jj-1))^2;
        dfiddi2FDFD(jj-1)=(-dfiddi(jj-1)+dfiddi(jj+1))/(di(jj)-di(jj-1));
    end
    figure
    subplot(3,1,1)
    plot(di,Fi)
    subplot(3,1,2)
    plot(di,dfiddi,di,dfiddiFD)
    subplot(3,1,3)
    plot(di,d2fiddi2,'r-',di(2:end-1),dfiddi2FDFD,'b-')
    
    
end

%% Surface Generation functions

function [h]=SurfFor2ndDeriv(X,Y,Z)

    c=(log(abs(Z)));
    h=surf(X,Y,Z,c,'linestyle','none');
    colormean=mean(c(isfinite(c)));
    colorstd=std(c(isfinite(c)));
    colormean
    colorstd
    caxis([colormean-colorstd,colormean+colorstd]);
    
end

%% Well Fuck it
function [derivtenscalc2]=ExtractDataForDerivatives(snaxel,snakposition,snakPosIndex,smearLengthEps,smearType)
    
    
    for ii=length(snakposition):-1:1
        neighSub=FindObjNum([],[snaxel(ii).snaxprec],snakPosIndex);
        
        derivtenscalc(ii).index=snakposition(ii).index;
        derivtenscalc(ii).snaxprec=snaxel(ii).snaxprec;
        derivtenscalc(ii).precsub=neighSub;
        % extracting data from preexisting arrays
        derivtenscalc(ii).Dg_i=snakposition(ii).vectornotnorm;
        derivtenscalc(ii).Dg_m=snakposition(neighSub).vectornotnorm;
        derivtenscalc(ii).p_i=snakposition(ii).coord;
        derivtenscalc(ii).p_m=snakposition(neighSub).coord;
        
        derivtenscalc(ii).g1_i=snakposition(ii).vertInit;
        derivtenscalc(ii).g1_m=snakposition(neighSub).vertInit;
    end
    switch smearType
        case 'length'
            for ii=length(snakposition):-1:1
                derivtenscalc(ii).d_i=snaxel(ii).d;
                derivtenscalc(ii).d_m=snaxel(neighSub).d;
                % calculating data
                % 1st type of smearing - General length smearing
                derivtenscalc(ii).normFi=sqrt(smearLengthEps^2+sum((derivtenscalc(ii).p_i-derivtenscalc(ii).p_m).^2));
                
            end
        case 'd'
            for ii=length(snakposition):-1:1
                derivtenscalc(ii).d_i=(1-2*smearLengthEps)*snaxel(ii).d+smearLengthEps;
                derivtenscalc(ii).d_m=(1-2*smearLengthEps)*snaxel(neighSub).d+smearLengthEps;
                % calculating data
                % 2nd type of smearing - Distance smearing
                derivtenscalc(ii).normFi=sqrt(sum(((derivtenscalc(ii).g1_i+derivtenscalc(ii).Dg_i*derivtenscalc(ii).d_i)-(derivtenscalc(ii).g1_m+derivtenscalc(ii).Dg_m*derivtenscalc(ii).d_m)).^2));
                
            end
    end
    
    
    for ii=length(snakposition):-1:1
        [derivtenscalc(ii).a_i,...
            derivtenscalc(ii).a_m,...
            derivtenscalc(ii).a_im,...
            derivtenscalc(ii).b_i,...
            derivtenscalc(ii).b_m,...
            derivtenscalc(ii).c]=...
            Calc_LengthDerivCoeff(...
            derivtenscalc(ii).Dg_i,derivtenscalc(ii).Dg_m,...
            derivtenscalc(ii).g1_i,derivtenscalc(ii).g1_m);
        
        [derivtenscalc2(ii)]=CalculateDerivatives(derivtenscalc(ii));
        
    end
    testnan=find(isnan([derivtenscalc2(:).d2fiddi2]));
    if ~isempty(testnan)
        testnan
    end
end

function [derivtenscalcII]=CalculateDerivatives(derivtenscalcII)
    
    varExtract={'a_i','a_m','a_im','b_i','b_m','c','normFi','d_i','d_m'};
    
    for ii=1:length(varExtract)
        eval([varExtract{ii},'=derivtenscalcII.(varExtract{ii});'])
    end
    
    [derivtenscalcII.dfiddi]=Calc_DFiDdi(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.dfiddm]=Calc_DFiDdm(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.d2fiddi2]=Calc_D2FiDdi2(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.d2fiddm2]=Calc_DFiDdm2(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    [derivtenscalcII.d2fiddim]=Calc_D2FiDdim(a_i,a_m,a_im,b_i,b_m,c,normFi,d_i,d_m);
    
end

function [dfiddi]=Calc_DFiDdi(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    
    dfiddi=(2*a_i*di+a_im*dm+b_i)/(2*(normFi));
    
    
end
function [dfiddm]=Calc_DFiDdm(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    dfiddm=(2*a_m*dm+a_im*di+b_m)/(2*(normFi));
end
function [d2fiddi2]=Calc_D2FiDdi2(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    
    d2fiddi2=((2*a_i*2*(normFi)^2)...
        -(2*a_i*di+a_im*dm+b_i)^2) ...
    /(4*((normFi)^3));

end
function [d2fiddm2]=Calc_DFiDdm2(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)

    d2fiddm2=((2*a_m*2*(normFi)^2)...
        -(2*a_m*dm+a_im*di+b_m)^2) ...
    /(4*(normFi)^3);  
end
function [d2fiddim]=Calc_D2FiDdim(a_i,a_m,a_im,b_i,b_m,c,normFi,di,dm)
    
    
    d2fiddim=((a_im*2*(normFi)^2)...
        -((2*a_m*dm+a_im*di+b_m)*(2*a_i*di+a_im*dm+b_i))) ...
        /(4*(normFi)^3);  
    
    
end



