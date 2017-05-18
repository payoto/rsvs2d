

%% Generate Solutions

syms a c1 c2 m l x 
syms f1 f2 f3 y
f1=-sqrt(1-(m*a-c1)^2)/m+c2;
f2=-sqrt(1-(m*a+c1)^2)/m+c2;
% f1=1-(m*a-c1)^2-m^2*c2^2;
% f2=1-(m*a+c1)^2-m^2*c2^2;
f3=2*a*c2-1/(2*m^2)*((m*a-c1)*sqrt(1-(c1-a*m)^2)-asin(c1-m*a)+(m*a+c1)*sqrt(1-(c1+a*m)^2)+asin(c1+m*a))-l;
y=sqrt(1-(m*x-c1)^2)/m+c2;
solSys=solve([f1 f2 f3],'returnconditions',true);
assume(l,'real');
assume(c2,'real');
assume(c1,'real');
assume(m,'real');
assume(a,'real');
assume(l~=0);
assume(a~=0);
solSys2=solve([f1 f2 f3],'returnconditions',true);


%% Output Latex code

[str]=MakeLatexEquationCode(y,'y(x)=','');
[str]=[str;MakeLatexEquationCode(f1,'\mathrm{end\ condition\ 1:\ }','=0')];
[str]=[str;MakeLatexEquationCode(f2,'\mathrm{end\ condition\ 2:\ }','=0')];
[str]=[str;MakeLatexEquationCode(f3,'\mathrm{Integral\ condition\ :\ }','=0')];


%sol No assumption


%sol assumption



[str2]=[MakeLatexEquationCode(solSys2.c1,'c1=','')];
[str2]=[str2;MakeLatexEquationCode(solSys2.l,'l=','')];
[str2]=[str2;MakeLatexEquationCode(solSys2.m,'m=','')];
[str2]=[str2;MakeLatexEquationCode(solSys2.conditions,'\mathrm{Condition\ :\ }','')];

%% Test conditions

syms cond1(a,c2) cond2(a,c2)

cond1(a,c2)=solSys2.conditions(1);
cond2(a,c2)=solSys2.conditions(2);

% cond1R=@(a,c2) (a^3*(2*asin(((4*a^2 + c2^4)^(1/2) - c2^2)/(2*a)) + ((1 - ((4*a^2 + c2^4)^(1/2) - c2^2)^2/(4*a^2))^(1/2)*((4*a^2 + c2^4)^(1/2) - c2^2))/a) + c2*((4*a^2 + c2^4)^(1/2) - c2^2)^2 ~= 0 | (4*a^2 + c2^4)^(1/2) == c2^2) & in(((4*a^2 + c2^4)^(1/2) - c2^2)/a^2, 'real');
% cond2R=@(a,c2) in(((4*a^2 + c2^4)^(1/2) + c2^2)/a^2, 'real') & ((4*a^2 + c2^4)^(1/2) + c2^2 == 0 | a^3*(2*asin(((4*a^2 + c2^4)^(1/2) + c2^2)/(2*a)) + ((1 - ((4*a^2 + c2^4)^(1/2) + c2^2)^2/(4*a^2))^(1/2)*((4*a^2 + c2^4)^(1/2) + c2^2))/a) ~= c2*((4*a^2 + c2^4)^(1/2) + c2^2)^2);

cond1R=@(a,c2)  0 < a.^2 + c2.^2 & 4.*a.*c2 ~= (2.*asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + 2.*a.*(1./(a.^2 + c2.^2)).^(1./2).*(1 - a.^2./(a.^2 + c2.^2)).^(1./2)).*(a.^2 + c2.^2);
cond2R=@(a,c2)  4.*a.*c2 + (2.*asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + 2.*a.*(1./(a.^2 + c2.^2)).^(1./2).*(1 - a.^2./(a.^2 + c2.^2)).^(1./2)).*(a.^2 + c2.^2) ~= 0 & 0 < a.^2 + c2.^2;


constr1=@(a,c2)  2.*a.*c2 - (a.^2 + c2.^2).*(asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + (a.*abs(c2))./(a.^2 + c2.^2));
constr2=@(a,c2)  2.*a.*c2 + (a.^2 + c2.^2).*(asin(a.*(1./(a.^2 + c2.^2)).^(1./2)) + (a.*abs(c2))./(a.^2 + c2.^2));


xa=linspace(0.01,1,101);
yc2=linspace(-20,5,101);

[Xa,Yc2]=meshgrid(xa,yc2);

iscond1=(cond1R(Xa,Yc2));
iscond2=(cond2R(Xa,Yc2));
constrval1=constr1(Xa,Yc2);
constrval2=constr2(Xa,Yc2);

figure
subplot(1,2,1)
surf(Xa,Yc2,double(iscond2))
view(0,90)
subplot(1,2,2)
s(1)=surf(Xa,Yc2,real(constrval2));
hold on
s(2)=surf(Xa,Yc2,real(constrval2)*0);
s(2).FaceColor=[1 0 0 ];
s(2).LineStyle='none';

%% 








