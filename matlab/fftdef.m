
%%%  fftdef.m --- defining the grids for the fft/split method
%%%
%%%  2000-2001, Thomas Busch    (thbusch@phys.ucc.ie)
%%%  Version 1.0
%%%
%%%  Change history:
%%%  14 Nov 01 Thomas Busch    - starting better documentation 
%%%
%%%  Notes: defining various grids for the split-operator/fft
%%%         method for solving the GPE
%%%
%%%         [X,DX,PX,DPX]=FFTDEF(XMAX,NGX) 
%%%         where XMAX is the maximum value for the X-coordinate 
%%%         and NGX the amount of points in that direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,dx,px,dpx]=fftdef(xmax,Ngx);

 % maximum values
 pxmax=pi*Ngx/(2*xmax);
 % spacing in position and momentum space
 dx=2*xmax/Ngx;   dpx=2*pxmax/Ngx;
 % grid vectors, position and momentum space
 x=(1:Ngx)'*dx-xmax;   pxn=((1:Ngx)*dpx-pxmax)';   
 % reordination needed for the fourier transform
 px(Ngx/2+2:Ngx,1)=pxn(1:Ngx/2-1,1);   px(1:Ngx/2+1,1)=pxn(Ngx/2:Ngx,1);