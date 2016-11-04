function [x_fin,F,LOG]=SUBopt3D_SNopt(x_curr,obj_fun,I,O)
Fglo=[]; xglo=[]; stopopt=false; SOL=[];
funcalls=0; gradcalls=0;
maj_iter=0; min_iter=0;
x_hist=[]; dir_hist=[];
snprint 'print.out'
snset('Nonderivative linesearch')
snsetr('Function precision',O.fun_prec)
snseti('Verify Level',O.grad_verify_level)
snseti('Major iterations limit',O.SNopt_max_iter)
snsetr('Major optimality tolerance',O.Optimality_tol)
snscreen on
snsummary 'summary.out'
xlow=-inf*ones(size(x_curr));
xupp=inf*ones(size(x_curr));
xmul=-ones(size(x_curr));
xstate=-ones(size(x_curr));

iGfun=reshape((1:length(1))'*(ones(1,length(x_curr))),[],1);
jGvar=reshape(ones(length(1),1)*(1:length(x_curr)),[],1);

snoptfun=@(x)snopt_fun_wrapper(x,obj_fun,O);

[x_fin,F,inform,xmul,Fmul] = snopt ( x_curr, xlow, xupp, xmul, xstate,...
    -inf, inf, 0, 1, snoptfun ,[],[],[],iGfun,jGvar);

if stopopt==false
    try
        Fglo(end+1,:)=F;
        xglo(:,end+1)=x_curr;
    catch
        F
        x_curr
    end
else
    F=Fglo(end,:)';
    x_curr=xglo(:,end)';
end
LOG.F=Fglo;
LOG.x=xglo;
LOG.SOL=SOL;
LOG.funcalls=funcalls;
LOG.gradcalls=gradcalls;
LOG.xhist=x_hist;
LOG.xmul=xmul;
LOG.Fmul=Fmul;
LOG.inform=inform;
snprint off
snsummary off
snend;

    function [F,G]=snopt_fun_wrapper(x,obj_fun,O)
        if stopopt==true
            F=Fglo(end,:)';
            G=zeros(length(x),size(F,2))';
        else
            if isempty(x_hist)
                min_dist=0;
                ind=[];
            else
                [ind,min_dist]=closest_x(x',x_hist);
            end
            
            if nargout>1 && min_dist==0
                MAJOR=1;
            else
                MAJOR=0;
            end
            
            if nargout>1
                [F,G,SOLi] = obj_fun(x,1);
                funcalls=funcalls+1;
                gradcalls=gradcalls+1;
            else
                [F,G,SOLi] = obj_fun(x,0);
                funcalls=funcalls+1;
            end
            
            if MAJOR
                maj_iter=maj_iter+1;
                min_iter=1;
                Fglo(end+1,1)=F;
                xglo(:,end+1)=x;
                switch O.conv_type
                    case 'prop_max'
                        slope=diff(log10(Fglo(:,1)),[],1);
                        if size(slope,1)>max(O.w,O.m)
                            avslope=zeros(size(slope));
                            for i=1:size(slope,2)
                                avslope(:,i)=smooth(slope(:,i),O.m);
                            end
                            avslope=avslope(O.m:end,:);
                            maxslope=max(abs(avslope));
                            if all(abs(mean(slope(end-O.w+1:end,:),1))<O.t*maxslope)
                                stopopt = true
                            end
                        end
                    case 'none'
                end

                
            else
                min_iter=min_iter+1;
                F_last=F;
            end
            
            
            
            if strcmpi(O.solver,'SU2')
                dir=['history/' num2str(I) '.' num2str(maj_iter) '.' num2str(min_iter)];
            dir_hist{end+1}=dir;
                if nargout>1
                    saveSU2Data(dir,O,F,x,G)
                else
                    saveSU2Data(dir,O,F,x)
                end
            end
            
            SOL{maj_iter}{min_iter}=SOLi;
            
            x_string=sprintf('% 1.16e ', x);
            fprintf(O.fid_xhist,'%s\n',x_string);
            x_hist(end+1,:)=x;
            F_string=[];
            for i=1:length(O.headers)
                F_string=[F_string sprintf([' %-' num2str(O.cw(i)) O.var_ftype(i)],...
                    eval(O.var{i}))];
            end
            F_string=strrep(F_string,'NaN','  -');
            fprintf(O.fid_Fhist,'%s\n',F_string);
            if MAJOR
                fprintf(O.fid_Fiter,'%s\n',F_string);
            fprintf(O.fid_xiter,'%s\n',x_string);
            G_string=sprintf('% 1.16e ', G);
            fprintf(O.fid_Giter,'%s\n',G_string);
            end
        end
        
    end
end