function pout= NonlinearFFT(u,L,alpha2)
%This function uses Fourier transformation space to solve a 
%PDE (Heat equation) with second order space derivative
%    / d^2(T)\
%    |--------|= -alpha^2 * kap^2 * uhat
%    \ dx ^2  /
%Which will be converted to a ODE
%Use following inputs for this function
%   rhsFDM(timestep size,
%           input function which will be transformed,
%           Half length of the domain in meter,
%           Diffusion coeffcient,
%           option)

Npx=length(u);
uhat=fft(u);
kap=(2*pi/L)*(-Npx/2:Npx/2-1);
kap=fftshift(kap);
duhat=-alpha2*(kap.^2).*uhat;
pout=duhat;
end