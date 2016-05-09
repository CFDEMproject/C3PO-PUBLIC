% postproc.m: loads the data necessary and saves the graphs 

clear;
clc;
close all;

x=linspace(1.5,6,10);
val=linspace(1,4.5,10);
sol=linspace(1,4.5,500);
var_=linspace(1,4.5,10);
varSol_=linspace(1,4.5,500);
r=linspace(1,6,500);

name=0.5;

for i=1:10
name=name+0.25;
file=sprintf("../c3po_dataStorage_particles/time_0_%2.2f_processor0_particleCenter.json",name)

dat=loadjson([file]);
field_data =getfield(dat,"fields_at_particle_centers" );

cell_data = struct2cell(field_data);
val(i)=cell_data{1};
var_(i)=cell_data{4};

i=i+1;
end

for i=1:500
 sol(i)= 1;
 varSol_(i)=1/(5*r(i)^3);
end


plot(x,val,'ob',r,sol,'-r')

set(gca,'Fontsize',14)
legend('CPPPO','Analytical solution')
xlabel('R_f/R')
ylabel('U_f/U')

uiwait

plot(x,var_,'ob',r,varSol_,'-r')
axis([1 6 0 0.05])
set(gca,'Fontsize',14)
legend('CPPPO','Analytical solution')
xlabel('R_f/R')
ylabel('<U^2>')

uiwait

varErr=0;
aveErr=0;
for i=1:10
 
 aveErr= abs(val(i)-1)/10 + aveErr;
 
 varErr= abs(var_(i)-1/(5*x(i)^3))/10 +varErr;
 
end


filename = "../errorPotential.txt";
fid = fopen (filename, "w");
str=sprintf("Favre average -> average error:  %f \n",aveErr);
fputs (fid, str);
str=sprintf("Favre variance -> average error:  %f \n",varErr);
fputs (fid, str);

fclose (fid);


clear;
clc;
close all;


