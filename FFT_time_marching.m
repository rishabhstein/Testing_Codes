function fpr=FFT_time_marching(beta,delta,rho,c,dt,diffuhat,fpr,f1,f2,f3,f4)



fpr=fft(fpr);
f1=fft(f1);
f2=fft(f2);
f3=fft(f3);
f4=fft(f4);

What=-f1+2*fpr+diffuhat*dt^2;

Hhat=(beta/(rho*c^2))*(2*fpr.^2-5*f1.^2+4*f2.^2-f3.^2);

Dhat=(delta*(dt^2)/((2*dt^3)*c^2))*(5*fpr-18*f1+24*f2-14*f3+3*f4);

pnext=What+Hhat+Dhat;

    fpr=pnext;

fpr=real(ifft(fpr));
end