function phinew=updatephi(phi,U,V,h,ghostnum,dt)
modification=false;
% 1 for Sussman;2 for WENO
method=2;
[Ny,Nx]=size(phi);
N=Nx-2*ghostnum;
M=Ny-2*ghostnum;
start=ghostnum+1;
endy=M+ghostnum;
endx=N+ghostnum;

%update in x direction
phistar=zeros(Ny,Nx);
for j=start:endy
    for i=start:endx

        if(method==1)
        
        if(U(j,i+1)<=0)

            phi_phalf=phi(j,i+1)-h/2*(1+U(j,i+1)*dt/h)*(phi(j,i+2)-phi(j,i))/h/2;

        else

            phi_phalf=phi(j,i)+h/2*(1-U(j,i+1)*dt/h)*(phi(j,i+1)-phi(j,i-1))/h/2;

        end

        
        if(U(j,i)<=0)

            phi_mhalf=phi(j,i)-h/2*(1+U(j,i)*dt/h)*(phi(j,i+1)-phi(j,i-1))/h/2;

        else
            phi_mhalf=phi(j,i-1)+h/2*(1-U(j,i)*dt/h)*(phi(j,i)-phi(j,i-2))/h/2;

        end

        end

        if(method==2)
            eps=0;
            d0=1/10;d1=3/5;d2=3/10;
            if(U(j,i+1)>=0)

                beta0=13/12*(phi(j,i-2)-2*phi(j,i-1)+phi(j,i))^2+1/4*(phi(j,i-2)-4*phi(j,i-1)+3*phi(j,i))^2;
                beta1=13/12*(phi(j,i-1)-2*phi(j,i)+phi(j,i+1))^2+1/4*(phi(j,i-1)-phi(j,i+1))^2;
                beta2=13/12*(phi(j,i)-2*phi(j,i+1)+phi(j,i+2))^2+1/4*(3*phi(j,i)-4*phi(j,i+1)+phi(j,i+2))^2;
                f0=1/3*phi(j,i-2)-7/6*phi(j,i-1)+11/6*phi(j,i);
                f1=-1/6*phi(j,i-1)+5/6*phi(j,i)+1/3*phi(j,i+1);
                f2=1/3*phi(j,i)+5/6*phi(j,i+1)-1/6*phi(j,i+2);

            else
                beta0=13/12*(phi(j,i+3)-2*phi(j,i+2)+phi(j,i+1))^2+1/4*(phi(j,i+3)-4*phi(j,i+2)+3*phi(j,i+1))^2;
                beta1=13/12*(phi(j,i+2)-2*phi(j,i+1)+phi(j,i))^2+1/4*(phi(j,i+2)-phi(j,i))^2;
                beta2=13/12*(phi(j,i+1)-2*phi(j,i)+phi(j,i-1))^2+1/4*(3*phi(j,i+1)-4*phi(j,i)+phi(j,i-1))^2;
                f0=1/3*phi(j,i+3)-7/6*phi(j,i+2)+11/6*phi(j,i+1);
                f1=-1/6*phi(j,i+2)+5/6*phi(j,i+1)+1/3*phi(j,i);
                f2=1/3*phi(j,i+1)+5/6*phi(j,i)-1/6*phi(j,i-1);
         
            end
            
            alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
            sum_alpha=alpha0+alpha1+alpha2;
            omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
            phi_phalf=omega0*f0+omega1*f1+omega2*f2;

            if(U(j,i)>=0)

                beta0=13/12*(phi(j,i-3)-2*phi(j,i-2)+phi(j,i-1))^2+1/4*(phi(j,i-3)-4*phi(j,i-2)+3*phi(j,i-1))^2;
                beta1=13/12*(phi(j,i-2)-2*phi(j,i-1)+phi(j,i))^2+1/4*(phi(j,i-2)-phi(j,i))^2;
                beta2=13/12*(phi(j,i-1)-2*phi(j,i)+phi(j,i+1))^2+1/4*(3*phi(j,i-1)-4*phi(j,i)+phi(j,i+1))^2;
                f0=1/3*phi(j,i-3)-7/6*phi(j,i-2)+11/6*phi(j,i-1);
                f1=-1/6*phi(j,i-2)+5/6*phi(j,i-1)+1/3*phi(j,i);
                f2=1/3*phi(j,i-1)+5/6*phi(j,i)-1/6*phi(j,i+1);

            else
                beta0=13/12*(phi(j,i+2)-2*phi(j,i+1)+phi(j,i))^2+1/4*(phi(j,i+2)-4*phi(j,i+1)+3*phi(j,i))^2;
                beta1=13/12*(phi(j,i+1)-2*phi(j,i)+phi(j,i-1))^2+1/4*(phi(j,i+1)-phi(j,i-1))^2;
                beta2=13/12*(phi(j,i)-2*phi(j,i-1)+phi(j,i-2))^2+1/4*(3*phi(j,i)-4*phi(j,i-1)+phi(j,i-2))^2;
                f0=1/3*phi(j,i+2)-7/6*phi(j,i+1)+11/6*phi(j,i);
                f1=-1/6*phi(j,i+1)+5/6*phi(j,i)+1/3*phi(j,i-1);
                f2=1/3*phi(j,i)+5/6*phi(j,i-1)-1/6*phi(j,i-2);
         
            end
            alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
            sum_alpha=alpha0+alpha1+alpha2;
            omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
            phi_mhalf=omega0*f0+omega1*f1+omega2*f2;

        end

        term1=(phi_phalf*U(j,i+1)-phi_mhalf*U(j,i))/h;
        term2=(U(j,i+1)-U(j,i))/h*phi(j,i);
        phistar(j,i)=phi(j,i)+dt*(-term1+term2);

    end
end

phinew=zeros(Ny,Nx);
%update in y direction
for j=start:endy
    for i=start:endx

        if(method==1)

            if(V(j+1,i)<=0)
                phistar_phalf=phistar(j+1,i)-h/2*(1+V(j+1,i)*dt/h)*(phistar(j+2,i)-phistar(j,i))/h/2;
            else
                phistar_phalf=phistar(j,i)+h/2*(1-V(j+1,i)*dt/h)*(phistar(j+1,i)-phistar(j-1,i))/h/2;
            end


            if(V(j,i)<=0)
                phistar_mhalf=phistar(j,i)-h/2*(1+V(j,i)*dt/h)*(phistar(j+1,i)-phistar(j-1,i))/h/2;
            else
                phistar_mhalf=phistar(j-1,i)+h/2*(1-V(j,i)*dt/h)*(phistar(j,i)-phistar(j-2,i))/h/2;
            end

        end

        if(method==2)
            eps=0;
            d0=1/10;d1=3/5;d2=3/10;
            if(V(j,i)>=0)

                beta0=13/12*(phi(j-2,i)-2*phi(j-1,i)+phi(j,i))^2+1/4*(phi(j-2,i)-4*phi(j-1,i)+3*phi(j,i))^2;
                beta1=13/12*(phi(j-1,i)-2*phi(j,i)+phi(j+1,i))^2+1/4*(phi(j-1,i)-phi(j+1,i))^2;
                beta2=13/12*(phi(j,i)-2*phi(j+1,i)+phi(j+2,i))^2+1/4*(3*phi(j,i)-4*phi(j+1,i)+phi(j+2,i))^2;
                f0=1/3*phi(j-2,i)-7/6*phi(j-1,i)+11/6*phi(j,i);
                f1=-1/6*phi(j-1,i)+5/6*phi(j,i)+1/3*phi(j+1,i);
                f2=1/3*phi(j,i)+5/6*phi(j+1,i)-1/6*phi(j+2,i);

            else
                beta0=13/12*(phi(j+3,i)-2*phi(j+2,i)+phi(j+1,i))^2+1/4*(phi(j+3,i)-4*phi(j+2,i)+3*phi(j+1,i))^2;
                beta1=13/12*(phi(j+2,i)-2*phi(j+1,i)+phi(j,i))^2+1/4*(phi(j+2,i)-phi(j,i))^2;
                beta2=13/12*(phi(j+1,i)-2*phi(j,i)+phi(j-1,i))^2+1/4*(3*phi(j+1,i)-4*phi(j,i)+phi(j-1,i))^2;
                f0=1/3*phi(j+3,i)-7/6*phi(j+2,i)+11/6*phi(j+1,i);
                f1=-1/6*phi(j+2,i)+5/6*phi(j+1,i)+1/3*phi(j,i);
                f2=1/3*phi(j+1,i)+5/6*phi(j,i)-1/6*phi(j-1,i);
         
            end
            
            alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
            sum_alpha=alpha0+alpha1+alpha2;
            omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
            phistar_phalf=omega0*f0+omega1*f1+omega2*f2;

            if(V(j,i+1)>=0)

                beta0=13/12*(phi(j-3,i)-2*phi(j-2,i)+phi(j-1,i))^2+1/4*(phi(j-3,i)-4*phi(j-2,i)+3*phi(j-1,i))^2;
                beta1=13/12*(phi(j-2,i)-2*phi(j-1,i)+phi(j,i))^2+1/4*(phi(j-2,i)-phi(j,i))^2;
                beta2=13/12*(phi(j-1,i)-2*phi(j,i)+phi(j+1,i))^2+1/4*(3*phi(j-1,i)-4*phi(j,i)+phi(j+1,i))^2;
                f0=1/3*phi(j-3,i)-7/6*phi(j-2,i)+11/6*phi(j-1,i);
                f1=-1/6*phi(j-2,i)+5/6*phi(j-1,i)+1/3*phi(j,i);
                f2=1/3*phi(j-1,i)+5/6*phi(j,i)-1/6*phi(j+1,i);

            else
                beta0=13/12*(phi(j+2,i)-2*phi(j+1,i)+phi(j,i))^2+1/4*(phi(j+2,i)-4*phi(j+1,i)+3*phi(j,i))^2;
                beta1=13/12*(phi(j+1,i)-2*phi(j,i)+phi(j-1,i))^2+1/4*(phi(j+1,i)-phi(j-1,i))^2;
                beta2=13/12*(phi(j,i)-2*phi(j-1,i)+phi(j-2,i))^2+1/4*(3*phi(j,i)-4*phi(j-1,i)+phi(j-2,i))^2;
                f0=1/3*phi(j+2,i)-7/6*phi(j+1,i)+11/6*phi(j,i);
                f1=-1/6*phi(j+1,i)+5/6*phi(j,i)+1/3*phi(j-1,i);
                f2=1/3*phi(j,i)+5/6*phi(j-1,i)-1/6*phi(j-2,i);
         
            end
            alpha0=d0/(beta0+eps)^2;alpha1=d1/(beta1+eps)^2;alpha2=d2/(beta2+eps)^2;
            sum_alpha=alpha0+alpha1+alpha2;
            omega0=alpha0/sum_alpha;omega1=alpha1/sum_alpha;omega2=alpha2/sum_alpha;
            phistar_mhalf=omega0*f0+omega1*f1+omega2*f2;

        end

        termstar1=(phistar_phalf*V(j+1,i)-phistar_mhalf*V(j,i))/h;
        termstar2=(V(j+1,i)-V(j,i))/h*phi(j,i);
        phinew(j,i)=phistar(j,i)+dt*(-termstar1+termstar2);
          
    end
end

%set the boundary condition(non-slip)
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

%level set iteration to make sure the absolute value of gradient level set be zero
if (modification==true)
phi_mod_new=phinew;
epsilon=1e-10;
maxstep=1;
%伪时间步
% dtau=0.8*h;
% S=@(x) x./sqrt(x.^2+epsilon);
% for t=1:maxstep 
%     phi_mod_old=phi_mod_new;  
%     for j=start:endy
%         for i=start:endx
%             phix=(phi_mod_old(j,i+1)-phi_mod_old(j,i-1))/h/2;
%             phiy=(phi_mod_old(j+1,i)-phi_mod_old(j-1,i))/h/2;
%             norm=abs(sqrt(phix^2+phiy^2));
%             phi_mod_new(j,i)=phi_mod_old(j,i)+dtau*S(norm)*(1-norm);
%         end
%     end  
% end
% phinew=phi_mod_new;
% end

end
