% postproc.m: loads the data necessary and saves the graphs 

clear;
clc;
close all;

x=linspace(0,1,100);
f=linspace(0,1,100);
F=linspace(0,1,100);

n=1
for i=1:100

 F(i) = 1/(8*n*n);

endfor

F(1)
%{
plot(x,F,'-b')
set(gca,'Fontsize',14)
legend('f','F')
xlabel('x')
ylabel(' ')

uiwait
%}
