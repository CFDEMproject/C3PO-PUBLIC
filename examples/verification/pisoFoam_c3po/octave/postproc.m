% postproc.m: loads the data necessary and saves the graphs 

clear;
clc;
close all;

x=linspace(1.5,6,10);
val1=linspace(1,4.5,10);
val2=linspace(1,4.5,10);
sol=linspace(1,4.5,500);
var1_=linspace(1,4.5,10);
var2_=linspace(1,4.5,10);
varSol_=linspace(1,4.5,500);
r=linspace(1.1,6,500);

name=0.5;

for i=1:10
 name=name+0.25;
 file=sprintf("../c3po_dataStorage_particles/time_0.25_%2.2f_processor1_particleCenter.h5",name)
 
 load (file);
 val1(i)= U_Filter_0(1,1);
 var1_(i)=U_Filter_variance0_0(1,1);
 close all ;
 
 file=sprintf("../c3po_dataStorage_particles/time_0.5_%2.2f_processor1_particleCenter.h5",name)
 load (file) ;
 val2(i)= U_Filter_0(1,1);
 var2_(i)=U_Filter_variance0_0(1,1);
 close all ;
end


for i=1:500
 sol(i)= 1 - 3/2 * (r(i)*r(i)-1)/(r(i)*r(i)*r(i)-1);
 varSol_(i)=(18*r(i)^5 -32*r(i)^4 + 14*r(i)^3-3*r(i)^2+2*r(i)+1)/(20*r(i)^7 + 40*r(i)^6 + 60*r(i)^5+40*r(i)^4+20*r(i)^3);
end


plot(x,val1,'ob',x,val2,'*g',r,sol,'-r')

set(gca,'Fontsize',14)
legend('CPPPO time=0.25','CPPPO time=0.5','Analytical solution',"southeast")
xlabel('R_f/R')
ylabel('U_f/U')

print('-depsc','verificationAve.eps')

plot(x,var1_,'ob',x,var2_,'*g',r,varSol_,'-r')

set(gca,'Fontsize',14)
legend('CPPPO time=0.25','CPPPO time=0.5','Analytical solution',"southeast")
xlabel('R_f/R')
ylabel('<U^2>')

print('-depsc','verificationVar.eps')


varErr=0;
aveErr=0;
for i=1:10
 
 avee=1 - 3/2 * (x(i)*x(i)-1)/(x(i)*x(i)*x(i)-1);
 vare=(18*x(i)^5 -32*x(i)^4 + 14*x(i)^3-3*x(i)^2+2*x(i)+1)/(20*x(i)^7 + 40*x(i)^6 + 60*x(i)^5+40*x(i)^4+20*x(i)^3);
 
 aveErr= abs(val2(i)-avee)/10 + aveErr;
 
 varErr= abs(var2_(i)-vare)/10 +varErr;
 
end


filename = "../error.txt";
fid = fopen (filename, "w");
str=sprintf("Favre average -> average error:  %f \n",aveErr);
fputs (fid, str);
str=sprintf("Favre variance -> average error:  %f \n",varErr);
fputs (fid, str);

fclose (fid);



clear;
clc;
close all;

