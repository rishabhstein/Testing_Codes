%This code is developed for learning purpose to solve non-linear
%PDEs like 1-D Westervelt equation using Fourier transformation 
%and solving the algebric equations with FDM method.

clear all; close all; clc;

feq=10^5;%^5;               %frequence of input signal in Hz
c=1481.44;              %Speed of sound (m/sec)
xsh=0.52;               %Shock traveLimitlesslling distance
Td=6/feq;               %Source envelope delay
Tend=2*Td;              %Source end time 2*Td;
Tw=3/feq;               %Source envelope width(Variance)
L=20;%0.95*xsh+Tend*c;  %Length of domain

nx=100000;                       %no. of elements
dx=L/nx;
%dt=0.0001;%10^-6;%
dt=2*dx/(pi*c);                 %time step
x=-L/2:dx:L/2-dx;               %Meshing

p0=1*10^6;              %Input maximum pressure (Pascal)
omega=2*pi*feq;         %Angular frequency
delta=1*4.1*10^-6;      %Attenuation coefficient
beta=1*10;              %non-linearity coefficient
rho=999.6;              %density

%Gaussian input signal
p=0*x;
pfr=0*p; %p for next    step  (n+1)
ppr=p;   %p for current step  (n)
pp1=0*p; %p for previous step (n-1)
pp2=0*p; %p for previous step (n-2)
pp3=0*p; %p for previous step (n-3)
pp4=0*p; %p for previous step (n-4)

plotx=x;

for time=1:50000
    t=time*dt; 
    %------Boundary Condition---%
    ppr(1)=p0*sin(omega*(t-Td))*exp(-((t-Td)/(Tw/2))^2);
    ppr(nx)=0;   
    %---------------------------%
    duhat=NonlinearFFT(ppr,L,c);
    pfr=FFT_time_marching(beta,delta,rho,c,dt,duhat,ppr,pp1,pp2,pp3,pp4);
%with attenuation and non-linearity terms

    pp4=pp3;
    pp3=pp2;
    pp2=pp1;
    pp1=ppr;
    ppr=pfr;
    pFDM=pfr;
     
%    if mod(time,10)==0       
%    plot(plotx,(pFDM)/p0);  
%    str = sprintf('time %f' , t);
%    title(str);
%    pause(dt);
%    end

end
plot(plotx,(pFDM)/p0);
    str = sprintf('time %f' , t);
    title(str);