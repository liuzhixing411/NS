clc;
clear;
close all;

%% 定义变量

L=1;
W=1;
Initial_x=[0.25,0.75,0.25,0.75];
Initial_y=[0.25,0.25,0.75,0.75];
r=[0.1,0.1,0.1,0.1];
bubble_num=length(r);
h=0.01;
dt=0.005;
t=1;
N=L/h;
M=W/h;
T=t/dt+1;
tol=1e-6;

phi=zeros(M,N,T);
F=zeros(M,N,T);
U=zeros(M+1,N+1,T);
V=zeros(M+1,N+1,T);

X=linspace(0,L,N);
Y=linspace(0,W,M);
%% 写一个更新phi的函数
function phi_modified=updatephi(phi,h)
M=size(phi,1);
N=size(phi,2);
phinew=phi;
epsilon=1e-6;
maxstep=1;
%伪时间步
dtau=0.001;
S=@(x) x./sqrt(x.^2+epsilon);
for t=1:maxstep 
    phiold=phinew;  
    for j=1:M
        for i=1:N
              
            %先计算phi的梯度
            %x方向差分

            if(i==1)
                phix=(phiold(j,i+1)-phiold(j,N))/h/2;
            elseif(i==N)
                phix=(phiold(j,1)-phiold(j,i-1))/h/2;
            else
                phix=(phiold(j,i+1)-phiold(j,i-1))/h/2;
            end
            %y方向差分
            if(j==1)
                phiy=(phiold(j+1,i)-phiold(M,i))/h/2;
            elseif(j==M)
                phiy=(phiold(1,i)-phiold(j-1,i))/h/2;
            else
                phiy=(phiold(j+1,i)-phiold(j-1,i))/h/2;
            end
            %计算phi梯度的绝对值
            norm=abs(sqrt(phix^2+phiy^2));
            phinew(j,i)=phiold(j,i)+dtau*S(norm)*(1-norm);
        end
    end  
end
phi_modified=phinew;
end

%% 写一个验证归一化体积守恒的函数
function rate=conservation(F)
M=size(F,1);
N=size(F,2);
T=size(F,3);
V0=0;
V1=0;
for j=1:M
    for i=1:N
        V0=V0+(1-F(j,i,1));
        V1=V1+(1-F(j,i,T));
    end
end
rate=(V1-V0)/V0;

end

%% 写一个对F truncation的函数
function F_tr=truncation(F,tol)
if(F>1-tol)
    F_tr=1;
elseif(F<tol)

    F_tr=0;
else
    F_tr=F;
end
end

%% 定义一个直线绘制函数
function plotline(a,b,c,xmin,xmax,ymin,ymax)
%定义数组存储
xy=zeros(2,4);

%左边界
y_lf=(-a*xmin-c)/b;
if((y_lf>=ymin)&&(y_lf<ymax))
    xy(1,1)=xmin;
    xy(2,1)=y_lf;
end

%下边界
x_down=(-b*ymin-c)/a;
if((x_down>xmin)&&(x_down<=xmax))
    xy(1,2)=x_down;
    xy(2,2)=ymin;
end

%右边界
y_right=(-a*xmax-c)/b;
if((y_right>ymin)&&(y_right<=ymax))
    xy(1,3)=xmax;
    xy(2,3)=y_right;
end

%上边界
x_up=(-b*ymax-c)/a;
if((x_up>=xmin)&&(x_up<xmax))
    xy(1,4)=x_up;
    xy(2,4)=ymax;
end

flag=0;
for j=1:4

    if(flag==2)
        break;
    end

    if((xy(1,j)~=0)&&(flag==0))
        x_start=xy(1,j);
        y_start=xy(2,j);
        flag=flag+1;
        continue;
    end

    if((xy(1,j)~=0)&&(flag==1))
        x_end=xy(1,j);
        y_end=xy(2,j);
        flag=flag+1;
    end

end

if(flag==2)
plot([x_start,x_end],[y_start,y_end],'r');
end

end
%% 写一个找拟合平面的函数
function [a,b,c]=Emin(phi,j,i,H,h)
        [M,N]=size(phi);
        A=zeros(3,3);
        P=zeros(3,1);
        
        for m=1:3
           termp=0;

           for p=-1:1
               for q=-1:1
           pp=p;
           qq=q;

           %防越界的周期条件

           if((j+p)>M)
               pp=1-M;
           elseif((j+p)<1)
               pp=M-1;
           end

           if((i+q)>N)
               qq=1-N;
           elseif((i+q)<1)
               qq=N-1;
           end
      

           if(m==1)
               if((p==0)&&(q==0))
                   termp=termp+16*H(phi(j+pp,i+qq))*phi(j+pp,i+qq)*qq*h;
               else
                   termp=termp+H(phi(j+pp,i+qq))*phi(j+pp,i+qq)*qq*h;
               end
           elseif(m==2)
               if((p==0)&&(q==0))
                   termp=termp+16*H(phi(j+pp,i+qq))*phi(j+pp,i+qq)*pp*h;
               else
                   termp=termp+H(phi(j+pp,i+qq))*phi(j+pp,i+qq)*pp*h;
               end
           else
               if((p==0)&&(q==0))
                   termp=termp+16*H(phi(j+pp,i+qq))*phi(j+pp,i+qq);
               else
                   termp=termp+H(phi(j+pp,i+qq))*phi(j+pp,i+qq);
               end
           end

               end
           end

           P(m,1)=termp;


           for col=1:3
               
               term=0;

               for p=-1:1
                   for q=-1:1
                       pp=p;
                       qq=q;
                       %防越界的周期条件
                       if((j+p)>M)
                           pp=1-M;
                       elseif((j+p)<1)
                           pp=M-1;
                       end

                       if((i+q)>N)
                           qq=1-N;
                       elseif((i+q)<1)
                           qq=N-1;
                       end
           


                       if(m==1)
                           
                           if(col==1)
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(q*h)*(q*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(q*h)*(q*h);
                               end
                           
                           elseif(col==2)
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(p*h)*(q*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(p*h)*(q*h);
                               end
                           
                           else
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(q*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(q*h);
                               end

                           end
                       

                       elseif(m==2)

                           if(col==1)
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(q*h)*(p*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(q*h)*(p*h);
                               end
                           
                           elseif(col==2)
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(p*h)*(p*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(p*h)*(p*h);
                               end
                           
                           else
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(p*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(p*h);
                               end

                           end

                       else

                           if(col==1)
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(q*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(q*h);
                               end
                           
                           elseif(col==2)
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq))*(p*h);
                               else
                                   term=term+H(phi(j+pp,i+qq))*(p*h);
                               end
                           
                           else
                               if((p==0)&&(q==0))
                                   term=term+16*H(phi(j+pp,i+qq));
                               else
                                   term=term+H(phi(j+pp,i+qq));
                               end

                           end

                       end

                   end

               end

               A(m,col)=term;
               
           end

           
        end

        solution=A\P;
        a=solution(1,1);
        b=solution(2,1);
        c=solution(3,1);
        norm=sqrt(a^2+b^2);
        a=a/norm;
        b=b/norm;
        c=c/norm;
        

end
%% 写一个用F求c的函数
function c_solution=backward(F,j,i,a_solution,b_solution,c,h)              
                
                m1=min(abs(a_solution),abs(b_solution));
                m2=max(abs(a_solution),abs(b_solution));
                m=m1/(m1+m2);
                V1=m/2/(1-m);

                if((F(j,i)>0)&&(F(j,i)<=V1))
                    alpha=sqrt(2*m*(1-m)*F(j,i));
                elseif(F(j,i)<=1/2)
                    alpha=F(j,i)*(1-m)+m/2;
                elseif(F(j,i)<=(1-V1))
                    alpha=m/2+F(j,i)*(1-m);
                else
                    alpha=1-sqrt(2*m*(1-m)*(1-F(j,i)));
                end
                
                c_solution=sign(c)*abs(alpha*h-h/2)*(m1+m2);
end

%% 写一个resolution可调的积分函数

%积分边界填相对中心的平移值
function Area=Integration(a,b,c,xmin,xmax,ymin,ymax)
if xmax < xmin, t = xmin; xmin = xmax; xmax = t; end
if ymax < ymin, t = ymin; ymin = ymax; ymax = t; end
M=1000;
N=1000;
xx=linspace(xmin,xmax,N+1);
yy=linspace(ymin,ymax,M+1);
dx=(xmax-xmin)/N;
dy=(ymax-ymin)/M;
count=0;
for p=1:M
    for q=1:N
        if((a*(xx(q)+xx(q+1))/2+b*(yy(p)+yy(p+1))/2+c)>0)
            count=count+1;
        end
    end
end
Area=count*dx*dy;

end

%% 定义界面平滑函数
epsilon=sqrt(2)*h;

dirac = @(d) ((d >= -epsilon) & (d <= epsilon)) .* (0.5/epsilon) .* (1 + cos(pi*d/epsilon));

%% 写一个返回上下，或右左通量的函数
function [F_phalf,F_mhalf]=Flux(F,phi,U,j,i,h,tol,dirac,N,dt)
%判断U*F+项
if(U(j,i+1)>=0)
    if(F(j,i)>tol && F(j,i)<1-tol)
        [a_solution,b_solution,c]=Emin(phi,j,i,dirac,h);
        c_solution=backward(F,j,i,a_solution,b_solution,c,h);
        Area=Integration(a_solution,b_solution,c_solution,h/2-U(j,i+1)*dt,h/2,-h/2,h/2);
        F_phalf=Area/abs(U(j,i+1)*dt*h);
    else
        F_phalf=F(j,i);
    end

else
    if(i==N)
        if((F(j,1) <tol  || F(j,1) >1-tol))
            F_phalf=F(j,1);
        else
            [a_solution,b_solution,c]=Emin(phi,j,1,dirac,h);
            c_solution=backward(F,j,1,a_solution,b_solution,c,h);
            Area=Integration(a_solution,b_solution,c_solution,-h/2,-h/2-U(j,i+1)*dt,-h/2,h/2);
            F_phalf=Area/abs(U(j,i+1)*dt*h);

        end

    else
        if((F(j,i+1) <tol  || F(j,i+1) >1-tol))
            F_phalf=F(j,i+1);
        else
            [a_solution,b_solution,c]=Emin(phi,j,i+1,dirac,h);
            c_solution=backward(F,j,i+1,a_solution,b_solution,c,h);
            Area=Integration(a_solution,b_solution,c_solution,-h/2,-h/2-U(j,i+1)*dt,-h/2,h/2);
            F_phalf=Area/abs(U(j,i+1)*dt*h);

        end
    end
end

%判断U*F-项
if(U(j,i)<=0)
    if(F(j,i)>tol && F(j,i)<1-tol)
        [a_solution,b_solution,c]=Emin(phi,j,i,dirac,h);
        c_solution=backward(F,j,i,a_solution,b_solution,c,h);
        Area=Integration(a_solution,b_solution,c_solution,-h/2,-h/2-U(j,i)*dt,-h/2,h/2);
        F_mhalf=Area/abs(U(j,i)*dt*h);
    else
        F_mhalf=F(j,i);
    end
else
    if(i==1)

        if(F(j,N) <tol || F(j,N) > 1-tol)
            F_mhalf=F(j,N);
        else
            [a_solution,b_solution,c]=Emin(phi,j,N,dirac,h);
            c_solution=backward(F,j,N,a_solution,b_solution,c,h);
            Area=Integration(a_solution,b_solution,c_solution,h/2-U(j,i)*dt,h/2,-h/2,h/2);
            F_mhalf=Area/abs(U(j,i)*dt*h);
        end

    else
        if(F(j,i-1) <tol || F(j,i-1) > 1-tol)
            F_mhalf=F(j,i-1);
        else
            [a_solution,b_solution,c]=Emin(phi,j,i-1,dirac,h);

            c_solution=backward(F,j,i-1,a_solution,b_solution,c,h);
            Area=Integration(a_solution,b_solution,c_solution,h/2-U(j,i)*dt,h/2,-h/2,h/2);
            F_mhalf=Area/abs(U(j,i)*dt*h);
        end

    end

end
end

%% 初始化

position=zeros(M,N);
%初始化phi

for j=1:M
    for i=1:N
        x=i*h-h/2;
        y=j*h-h/2;
        for k=1:bubble_num
        middle(k)=sqrt((x-Initial_x(k))^2+(y-Initial_y(k))^2)-r(k);
        
        left(k)=sqrt((x-Initial_x(k)+L)^2+(y-Initial_y(k))^2)-r(k);
        
        right(k)=sqrt((x-Initial_x(k)-L)^2+(y-Initial_y(k))^2)-r(k);
        
        up(k)=sqrt((x-Initial_x(k))^2+(y-Initial_y(k)-W)^2)-r(k);
        
        down(k)=sqrt((x-Initial_x(k))^2+(y-Initial_y(k)+W)^2)-r(k);
        
        left_down(k)=sqrt((x-Initial_x(k)+L)^2+(y-Initial_y(k)+W)^2)-r(k);
       
        left_up(k)=sqrt((x-Initial_x(k)+L)^2+(y-Initial_y(k)-W)^2)-r(k);
        
        right_up(k)=sqrt((x-Initial_x(k)-L)^2+(y-Initial_y(k)-W)^2)-r(k);
        
        right_down(k)=sqrt((x-Initial_x(k)-L)^2+(y-Initial_y(k)+W)^2)-r(k);
        end
             
        A=[middle;left;right;up;down;left_up;left_down;right_down;right_up];
        m=min(A);
        [phi(j,i,1),position(j,i)]=min(m);
        
        
    end
end




%初始化F
for j=1:M
    for i=1:N
        x=i*h-h/2;
        y=j*h-h/2;
        pos=position(j,i);
        if(phi(j,i,1)>=(h/sqrt(2)))
            F(j,i,1)=1;
        elseif(phi(j,i,1)<=(-h/sqrt(2)))
            F(j,i,1)=0;
        else
            %进一步切分
            resM=200;
            resN=200;
            count=0;
            xx=linspace(x-h/2,x+h/2,resM+1);
            yy=linspace(y-h/2,y+h/2,resN+1);


            for l=1:resM
                for e=1:resN
                    if(sqrt(((xx(l)+xx(l+1))/2-Initial_x(pos))^2+((yy(e)+yy(e+1))/2-Initial_y(pos))^2)>r(pos))
                        count=count+1;
                    end
                end
            end


            F(j,i,1)=(count/(resM*resN));

        end
    end
end



%速度场定义
for n=1:T
for j=1:(M+1)
        for i=1:(N+1)
            U(j,i,n)=-cos((i-0.5)*h)*sin((j-0.5)*h)*exp(-2*(n-1)*dt) ;
        end
end

for j=1:(M+1)
        for i=1:(N+1)
            V(j,i,n)=sin((i-0.5)*h)*cos((j-0.5)*h)*exp(-2*(n-1)*dt);
        end
end

end


%% 可视化
figure


%% 主循环

for n=1:(T-1)
    %通过phi计算phi_star
    phi_star=zeros(M,N);
    for j=1:M
        for i=1:N
            %判断U*phi+项
            if(U(j,i+1,n)<=0)

                if(i==(N-1))
                    phi_phalf=phi(j,i+1,n)-h/2*(1+U(j,i+1,n)*dt/h)*(phi(j,1,n)-phi(j,i,n))/h/2;
                elseif(i==N)
                    phi_phalf=phi(j,1,n)-h/2*(1+U(j,i+1,n)*dt/h)*(phi(j,2,n)-phi(j,i,n))/h/2;
                else                    
                    phi_phalf=phi(j,i+1,n)-h/2*(1+U(j,i+1,n)*dt/h)*(phi(j,i+2,n)-phi(j,i,n))/h/2;
                end

            else

                if(i==N)
                    phi_phalf=phi(j,i,n)+h/2*(1-U(j,1,n)*dt/h)*(phi(j,1,n)-phi(j,i-1,n))/h/2;
                elseif(i==1)
                    phi_phalf=phi(j,i,n)+h/2*(1-U(j,i+1,n)*dt/h)*(phi(j,i+1,n)-phi(j,N,n))/h/2;
                else
                    phi_phalf=phi(j,i,n)+h/2*(1-U(j,i+1,n)*dt/h)*(phi(j,i+1,n)-phi(j,i-1,n))/h/2;
                end

            end

            %判断U*phi-项
            if(U(j,i,n)<=0)
                if(i==N)
                    phi_mhalf=phi(j,i,n)-h/2*(1+U(j,i,n)*dt/h)*(phi(j,1,n)-phi(j,i-1,n))/h/2;
                elseif(i==1)
                    phi_mhalf=phi(j,i,n)-h/2*(1+U(j,i,n)*dt/h)*(phi(j,i+1,n)-phi(j,N,n))/h/2;
                else
                phi_mhalf=phi(j,i,n)-h/2*(1+U(j,i,n)*dt/h)*(phi(j,i+1,n)-phi(j,i-1,n))/h/2;
                end
            else
                if(i==2)
                    phi_mhalf=phi(j,i-1,n)+h/2*(1-U(j,i,n)*dt/h)*(phi(j,i,n)-phi(j,N,n))/h/2;
                elseif(i==1)
                    phi_mhalf=phi(j,N,n)+h/2*(1-U(j,i,n)*dt/h)*(phi(j,i,n)-phi(j,N-1,n))/h/2;

                else
                phi_mhalf=phi(j,i-1,n)+h/2*(1-U(j,i,n)*dt/h)*(phi(j,i,n)-phi(j,i-2,n))/h/2;
                end
            end


            term1=(phi_phalf*U(j,i+1,n)-phi_mhalf*U(j,i,n))/h;
            term2=(U(j,i+1,n)-U(j,i,n))/h*phi(j,i,n);
            phi_star(j,i)=phi(j,i,n)+dt*(-term1+term2);
            

        end
    end


    %计算Fstar
    
    F_star=zeros(M,N);
    for j=1:M
        for i=1:N

            [F_phalf,F_mhalf]=Flux(F(:,:,n),phi(:,:,n),U(:,:,n),j,i,h,tol,dirac,N,dt);       
            term1=(F_phalf*U(j,i+1,n)-F_mhalf*U(j,i,n))/h;
            term2=(U(j,i+1,n)-U(j,i,n))/h*F(j,i,n);
            F_star(j,i)=truncation(F(j,i,n)+dt*(-term1+term2),tol);
            
        end
    end



    %phi_n+1
    for j=1:M
        for i=1:N
            %判断V*phi+项
            if(V(j+1,i,n)<=0)
                if(j==(M-1))
                    phistar_phalf=phi_star(j+1,i)-h/2*(1+V(j+1,i,n)*dt/h)*(phi_star(1,i)-phi_star(j,i))/h/2;

                elseif(j==M)
                    phistar_phalf=phi_star(1,i)-h/2*(1+V(j+1,i,n)*dt/h)*(phi_star(2,i)-phi_star(j,i))/h/2;
                else
                    phistar_phalf=phi_star(j+1,i)-h/2*(1+V(j+1,i,n)*dt/h)*(phi_star(j+2,i)-phi_star(j,i))/h/2;
                end

            else
                if(j==M)
                    phistar_phalf=phi_star(j,i)+h/2*(1-V(j+1,i,n)*dt/h)*(phi_star(1,i)-phi_star(j-1,i))/h/2;
                elseif(j==1)
                    phistar_phalf=phi_star(j,i)+h/2*(1-V(j+1,i,n)*dt/h)*(phi_star(j+1,i)-phi_star(M,i))/h/2;
                else
                    phistar_phalf=phi_star(j,i)+h/2*(1-V(j+1,i,n)*dt/h)*(phi_star(j+1,i)-phi_star(j-1,i))/h/2;
                end
            end

            %判断V*phi-项
            if(V(j,i,n)<=0)

                if(j==M)
                    phistar_mhalf=phi_star(j,i)-h/2*(1+V(j,i,n)*dt/h)*(phi_star(1,i)-phi_star(j-1,i))/h/2;
                elseif(j==1)
                    phistar_mhalf=phi_star(j,i)-h/2*(1+V(j,i,n)*dt/h)*(phi_star(j+1,i)-phi_star(M,i))/h/2;
                else
                    phistar_mhalf=phi_star(j,i)-h/2*(1+V(j,i,n)*dt/h)*(phi_star(j+1,i)-phi_star(j-1,i))/h/2;
                end

            else

                if(j==2)
                    phistar_mhalf=phi_star(j-1,i)+h/2*(1-V(j,i,n)*dt/h)*(phi_star(j,i)-phi_star(M,i))/h/2;
                elseif(j==1)
                    phistar_mhalf=phi_star(M,i)+h/2*(1-V(j,i,n)*dt/h)*(phi_star(j,i)-phi_star(M-1,i))/h/2;
                else
                    phistar_mhalf=phi_star(j-1,i)+h/2*(1-V(j,i,n)*dt/h)*(phi_star(j,i)-phi_star(j-2,i))/h/2;
                end

            end


            termstar1=(phistar_phalf*V(j+1,i,n)-phistar_mhalf*V(j,i,n))/h;
            termstar2=(V(j+1,i,n)-V(j,i,n))/h*phi(j,i,n);
            phi(j,i,n+1)=phi_star(j,i)+dt*(-termstar1+termstar2);

        end
    end


    %计算F_n+1
    for j=1:M
        for i=1:N
            [Fstar_phalf,Fstar_mhalf]=Flux(F_star',phi_star',V(:,:,n)',i,j,h,tol,dirac,M,dt);
            
            termstar1=(Fstar_phalf*V(j+1,i,n)-Fstar_mhalf*V(j,i,n))/h;
            termstar2=(V(j+1,i,n)-V(j,i,n))/h*F(j,i,n);
             F(j,i,n+1)=truncation(F_star(j,i)+dt*(-termstar1+termstar2),tol);

        end
    end

   
    %用F反向修正phi

    for j=1:M
        for i=1:N
            if((abs(phi(j,i,n+1))<h/sqrt(2))&&(F(j,i,n+1)>tol)&&(F(j,i,n+1)<1-tol))
                [a_solution,b_solution,c]=Emin(phi(:,:,n+1),j,i,dirac,h);
                c_solution=backward(F(:,:,n+1),j,i,a_solution,b_solution,c,h);
                phi(j,i,n+1)=sign(phi(j,i,n+1))*abs(c_solution);
            end
        end
    end


    %phi的梯度修正,每隔几个循环修正一次
    
    phi(:,:,n+1)=updatephi(phi(:,:,n+1),h);

   
    clf;
    hold on
    axis equal tight   
    contour(X,Y,phi(:,:,n),[0 0],'b'); 
    contour(X,Y,F(:,:,n),[0.5 0.5],'r');
    
    axis([0 1 0 1]); 
    drawnow;
    hold off
    
end

matrixF=F(:,:,31);
matrixphi=phi(:,:,31);
