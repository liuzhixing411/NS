%% NS solver
function [unew,vnew,p]=NS(u,v,h,ghostnum,phi,dt,rhol,rhog,mul,mug,sigma)
[Ny,Nx]=size(phi);
N=Nx-2*ghostnum;
M=Ny-2*ghostnum;
start=ghostnum+1;
endy=M+ghostnum;
endx=N+ghostnum;

upwindorder=2;

g=9.8;
%function definition
epsilon=sqrt(2)*h;
dirac = @(x) ((x >= -epsilon) & (x <= epsilon)) .* (0.5./epsilon) .* (1 + cos(pi.*x./epsilon));
H = @(x) (x < -epsilon).*0 + ((x >= -epsilon) & (x <= epsilon)).*0.5.*(1 + x./epsilon + sin(pi.*x./epsilon)./pi) + (x > epsilon).*1;
rho=@(x) rhol.*H(x)+rhog.*(1-H(x));
mu=@(x) mul.*H(x)+mug.*(1-H(x));

%% Step1 get ustar

% advection term in x direction
phihalfx=zeros(Ny+1,Nx+1);

for j=start:endy+1
    for i=start:endx+1
        phihalfx(j,i)=(phi(j,i)+phi(j,i-1))/2;
    end
end

adux=zeros(Ny+1,Nx+1);
aduy=zeros(Ny+1,Nx+1);

if (upwindorder==1)
for j=start:endy+1
    for i=start:endx+1
        if(u(j,i)>0)
            adux(j,i)=u(j,i)*(u(j,i)-u(j,i-1))/h;
        else
            adux(j,i)=u(j,i)*(u(j,i+1)-u(j,i))/h;
        end

        if(v(j,i)>0)
            aduy(j,i)=v(j,i)*(u(j,i)-u(j-1,i))/h;
        else
            aduy(j,i)=v(j,i)*(u(j+1,i)-u(j,i))/h;
        end
    end
end
end

if (upwindorder==2)
for j=start:endy+1
    for i=start:endx+1
        if(u(j,i)>0)
            adux(j,i)=u(j,i)*(3*u(j,i)-4*u(j,i-1)+u(j,i-2))/(2*h);
        else
            adux(j,i)=u(j,i)*(-3*u(j,i)+4*u(j,i+1)-u(j,i+2))/(2*h);
        end

        if(v(j,i)>0)
            aduy(j,i)=v(j,i)*(3*u(j,i)-4*u(j-1,i)+u(j-2,i))/(2*h);
        else
            aduy(j,i)=v(j,i)*(-3*u(j,i)+4*u(j+1,i)-u(j+2,i))/(2*h);
        end
    end
end
end

adu=adux+aduy;

%diffusion term in x direction
diffux=zeros(Ny+1,Nx+1);
diffuy=zeros(Ny+1,Nx+1);
diffvx=zeros(Ny+1,Nx+1);
diffvy=zeros(Ny+1,Nx+1);
diffu=zeros(Ny+1,Nx+1);
diffv=zeros(Ny+1,Nx+1);

for j=start:endy+1
    for i=start:endx+1

        if(u(j,i)>0)
            diffux(j,i)=(u(j,i)-u(j,i-1))/h;
            diffvx(j,i)=(v(j,i)-v(j,i-1))/h;
        else
            diffux(j,i)=(u(j,i+1)-u(j,i))/h;
            diffvx(j,i)=(v(j,i+1)-v(j,i))/h;
        end

        if(v(j,i)>0)
            diffuy(j,i)=(u(j,i)-u(j-1,i))/h;
            diffvy(j,i)=(v(j,i)-v(j-1,i))/h;
        else
            diffuy(j,i)=(u(j+1,i)-u(j,i))/h;
            diffvy(j,i)=(v(j+1,i)-v(j,i))/h;
        end
    end  
end

for j=start:endy+1
    for i=start:endx+1

        diff1=mu(phihalfx(j,i))*((u(j+1,i)-2*u(j,i)+u(j-1,i))/h^2+(u(j,i+1)-2*u(j,i)+u(j,i-1))/h^2);

        diff2x=((mu(phi(j,i))-mu(phi(j,i-1)))/h)*2*diffux(j,i);

        diffvx_mod=(diffvx(j,i)+diffvx(j,i+1)+diffvx(j+1,i)+diffvx(j+1,i+1))/4;

        diff2y=((mu(phi(j,i))-mu(phi(j-1,i)))/h)*(diffuy(j,i)+diffvx_mod);

        diffu(j,i)=(diff1+diff2x+diff2y)/(rho(phihalfx(j,i)));
    end
end

%surface tension term
nux=zeros(Ny+1,Nx+1);
nvy=zeros(Ny+1,Nx+1);
curu=zeros(Ny+1,Nx+1);
curv=zeros(Ny+1,Nx+1);
nx=zeros(Ny+1,Nx+1);
ny=zeros(Ny+1,Nx+1);
cur=zeros(Ny+1,Nx+1);
surfu=zeros(Ny+1,Nx+1);
surfv=zeros(Ny+1,Nx+1);

%normal vector
for j=start-1:endy+2
    for i=start-1:endx+2

        nxx=(phi(j,i+1)-phi(j,i-1))/h/2;
        nyy=(phi(j+1,i)-phi(j-1,i))/h/2;
        norm=sqrt(nxx^2+nyy^2);
        nx(j,i)=nxx/norm;
        ny(j,i)=nyy/norm;

    end
end

%curvature
for j=start:endy+1
    for i=start:endx+1
        gradnx=(nx(j,i+1)-nx(j,i-1))/h/2;
        gradny=(ny(j+1,i)-ny(j-1,i))/h/2;
        cur(j,i)=gradny+gradnx;
    end
end

for j=start:endy+1
    for i=start:endx+1
        curu(j,i)=(cur(j,i)*abs(phi(j,i-1))+cur(j,i-1)*abs(phi(j,i)))/(abs(phi(j,i))+abs(phi(j,i-1)));
        curv(j,i)=(cur(j,i)*abs(phi(j-1,i))+cur(j-1,i)*abs(phi(j,i)))/(abs(phi(j,i))+abs(phi(j-1,i)));
        nux(j,i)=(nx(j,i)+nx(j,i-1))/2;
        nvy(j,i)=(ny(j,i)+ny(j-1,i))/2;
    end
end

%solve surface tension
for j=start:endy+1
    for i=start:endx+1
        %surfu(j,i)=sigma*curu(j,i)*nux(j,i)*dirac(phihalfx(j,i))/(rho(phihalfx(j,i)));
        
        surfu(j,i)=sigma*curu(j,i)/(rho(phihalfx(j,i)))*(H(phi(j,i))-H(phi(j,i-1)))/h;

    end
end



%solve ustar
ustar=u+(-adu+diffu-surfu)*dt;

%set boundary condition(non-slip)
for j=1:Ny+1
    for i=1:ghostnum
        ustar(j,i)=-ustar(j,2*ghostnum+2-i);
    end
end

for j=1:Ny+1
    for i=endx+1:Nx+1
        ustar(j,i)=-ustar(j,2*endx+2-i);
    end
end

for j=1:ghostnum
    for i=1:Nx+1
        ustar(j,i)=-ustar(2*ghostnum+1-j,i);
    end
end

for j=endy+1:Ny+1
    for i=1:Nx+1
        ustar(j,i)=-ustar(2*endy+1-j,i);
    end
end

ustar(:,ghostnum+1)=zeros(size(ustar(:,ghostnum+1)));
ustar(:,endx+1)=zeros(size(ustar(:,endx+1)));

%% Step2 get vstar

% advection term in y direction
phihalfy=zeros(Ny+1,Nx+1);

for j=start:endy+1
    for i=start:endx+1
        phihalfy(j,i)=(phi(j,i)+phi(j-1,i))/2;
    end
end

advx=zeros(Ny+1,Nx+1);
advy=zeros(Ny+1,Nx+1);

if(upwindorder==1)
for j=start:endy+1
    for i=start:endx+1
        if(u(j,i)>0)
            advx(j,i)=u(j,i)*(v(j,i)-v(j,i-1))/h;
        else
            advx(j,i)=u(j,i)*(v(j,i+1)-v(j,i))/h;
        end

        if(v(j,i)>0)
            advy(j,i)=v(j,i)*(v(j,i)-v(j-1,i))/h;
        else
            advy(j,i)=v(j,i)*(v(j+1,i)-v(j,i))/h;
        end
    end
end

if (upwindorder==2)
for j=start:endy+1
    for i=start:endx+1
        if(u(j,i)>0)
            advx(j,i)=u(j,i)*(3*v(j,i)-4*v(j,i-1)+v(j,i-2))/(2*h);
        else
            advx(j,i)=u(j,i)*(-3*v(j,i)+4*v(j,i+1)-v(j,i+2))/(2*h);
        end

        if(v(j,i)>0)
            advy(j,i)=v(j,i)*(3*v(j,i)-4*v(j-1,i)+v(j-2,i))/(2*h);
        else
            advy(j,i)=v(j,i)*(-3*v(j,i)+4*v(j+1,i)-v(j+2,i))/(2*h);
        end
    end
end
end
end

adv=advx+advy;

%diffusion term in y direction
for j=start:endy+1
    for i=start:endx+1

        diff1=mu(phihalfy(j,i))*((v(j+1,i)-2*v(j,i)+v(j-1,i))/h^2+(v(j,i+1)-2*v(j,i)+v(j,i-1))/h^2);
        diff2x=((mu(phi(j,i))-mu(phi(j-1,i)))/h)*2*diffvy(j,i);

        diffuy_mod=(diffuy(j,i)+diffuy(j,i+1)+diffuy(j+1,i)+diffuy(j+1,i+1))/4;

        diff2y=((mu(phi(j,i))-mu(phi(j,i-1)))/h)*(diffuy_mod+diffvx(j,i));

        diffv(j,i)=(diff1+diff2x+diff2y)/(rho(phihalfy(j,i)));
    end
end

%solve surface tension term

for j=start:endy+1
    for i=start:endx+1
        %surfv(j,i)=sigma*curv(j,i)*nvy(j,i)*dirac(phihalfy(j,i))/(rho(phihalfy(j,i)));
        
        surfv(j,i)=sigma*curv(j,i)/(rho(phihalfy(j,i)))*(H(phi(j,i))-H(phi(j-1,i)))/h;
    end
end

%solve vstar
vstar=v+(-adv+diffv-surfv-g)*dt;

%set boundary condition(non-slip)
for j=1:Ny+1
    for i=1:ghostnum
        vstar(j,i)=-vstar(j,2*ghostnum+1-i);
    end
end

for j=1:Ny+1
    for i=endx+1:Nx+1
        vstar(j,i)=-vstar(j,2*endx+1-i);
    end
end

for j=1:ghostnum
    for i=1:Nx+1
        vstar(j,i)=-vstar(2*ghostnum+2-j,i);
    end
end

for j=endy+1:Ny+1
    for i=1:Nx+1
        vstar(j,i)=-vstar(2*endy+2-j,i);
    end
end

vstar(ghostnum+1,:)=zeros(size(vstar(ghostnum+1,:)));
vstar(endy+1,:)=zeros(size(vstar(endy+1,:)));

%Use Possion solver to modify velocity
[unew, vnew,p]=Possion(rho(phi), ustar, vstar, h, h, dt, ghostnum);

end
