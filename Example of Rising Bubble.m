clc;
clear;
close all
%% Parameters
L=0.1;
W=0.1;
h=0.0005;
timer=0;
x_c=0.05;
y_c=0.025;
r=0.01;
T=500;
M=L/h;
N=W/h;
ghostnum=3;
start=ghostnum+1;
endy=M+ghostnum;
endx=N+ghostnum;
Ny=ghostnum*2+M;
Nx=ghostnum*2+N;
phi=zeros(M,N);
rhol=1000;
rhog=500;
mul=1e-3;
mug=1.5e-5;
sigma=7.28e-2;
phiold=zeros(Ny,Nx);
phinew=zeros(Ny,Nx);

uold=zeros(Ny+1,Nx+1);
unew=zeros(Ny+1,Nx+1);
vold=zeros(Ny+1,Nx+1);
vnew=zeros(Ny+1,Nx+1);

X=linspace(0,W,N);
Y=linspace(0,L,M);


%% Initialization
for j=1:M
    for i=1:N
        x=i*h-h/2;
        y=j*h-h/2;       
        phi(j,i)=sqrt((x-x_c)^2+(y-y_c)^2)-r;
    end
end

phinew(start:endy,start:endx)=phi;

for j=1:ghostnum
    for i=1:Nx
        phinew(j,i)=phinew(2*ghostnum+1-j,i);
    end
end

for j=endy+1:Ny
    for i=1:Nx
        phinew(j,i)=phinew(2*endy+1-j,i);
    end
end

for j=1:Ny
    for i=1:ghostnum
        phinew(j,i)=phinew(j,2*ghostnum+1-i);
    end
end

for j=1:Ny
    for i=endx+1:Nx
        phinew(j,i)=phinew(j,2*endx+1-i);
    end
end

%% main loop

for n=1:T
    phiold=phinew;
    uold=unew;
    vold=vnew;
    %introduce dynamic time step based on CFL
    eps=0.01;
    Umax=max(max(uold));
    Vmax=max(max(vold));
    Max=max(Umax,Vmax);
    dt=min(0.25*h/Max,eps);
    timer=timer+dt;

    phinew=updatephi(phiold,uold,vold,h,ghostnum,dt);
    [unew,vnew,p]=NS(uold,vold,h,ghostnum,phiold,dt,rhol,rhog,mul,mug,sigma);
    hold on
    clf
    contour(X,Y,phinew(start:endy,start:endx),[0 0],'b');
    %contour(X,Y,vnew(start:endy,start:endx),40,'r');
    axis equal;
    pause(0.1);
    drawnow;
end



