
clear; close all; clc;

%% Parameters
Lx = 8; Ly = 8;  
Nx = 160; Ny = 160;   
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
h = x(2)-x(1);
[X, Y] = meshgrid(x, y);


xc1 = 4; yc1 = 4;    
r1 = 1;              

xc2=6;yc2=6;
r2=0.1;

% f(x) it's distance field, so it should be 1
f = ones(Ny, Nx);    

%% Initialization, unknown nodes need to be set as infinity

u = Inf(Ny, Nx);

% Use nodes on the interface as boundary condition (Dirichlet)
phi = min(sqrt((X - xc1).^2 + (Y - yc1).^2)-r1,sqrt((X - xc2).^2 + (Y - yc2).^2)-r2);
boundary_mask = abs(phi) <= (h*0.6);  
for j=1:Ny
    for i=1:Nx
        if(boundary_mask(j,i)==1)
            u(j,i)=phi(j,i);
        end
    end
end


%% Fast Sweeping 
%do 4 directions in turn
tol = 1e-6;
max_sweeps = 4;   % max sweep time, it's better to be 2 to 4 times
fprintf('Starting Fast Sweeping: Nx=%d Ny=%d, h=%g\n', Nx, Ny, h);

for sweep = 1:max_sweeps
    u_old = u;
    % 1) sweep: i=1->Nx, j=1->Ny
    for j = 1:Ny
        for i = 1:Nx
            if boundary_mask(j,i); continue; end
            if(i==1)
                a = min(u(j, Nx), u(j, i+1));
            elseif(i==Nx)
                a = min(u(j, i-1), u(j, 1));
            else
                a = min(u(j, i-1), u(j, i+1));
            end

            if(j==1)
                b = min(u(Ny, i), u(j+1, i));

            elseif(j==Ny)
                b = min(u(j-1, i), u(1, i));
            else
                b = min(u(j-1, i), u(j+1, i));
            end    % y-direction neighbors
            u(j,i) = local_update(a, b, f(j,i), h, u(j,i));
        end
    end
    % 2) sweep: i=Nx->1, j=1->Ny
    for j = 1:Ny
        for ii = Nx:-1:1
            i = ii;
            if boundary_mask(j,i); continue; end
            if(i==1)
                a = min(u(j, Nx), u(j, i+1));
            elseif(i==Nx)
                a = min(u(j, i-1), u(j, 1));
            else
                a = min(u(j, i-1), u(j, i+1));
            end

            if(j==1)
                b = min(u(Ny, i), u(j+1, i));

            elseif(j==Ny)
                b = min(u(j-1, i), u(1, i));
            else
                b = min(u(j-1, i), u(j+1, i));
            end
            u(j,i) = local_update(a, b, f(j,i), h, u(j,i));
        end
    end
    % 3) sweep: i=1->Nx, j=Ny->1
    for jj = Ny:-1:1
        for i = 1:Nx
            j = jj;
            if boundary_mask(j,i); continue; end
            if(i==1)
                a = min(u(j, Nx), u(j, i+1));
            elseif(i==Nx)
                a = min(u(j, i-1), u(j, 1));
            else
                a = min(u(j, i-1), u(j, i+1));
            end

            if(j==1)
                b = min(u(Ny, i), u(j+1, i));

            elseif(j==Ny)
                b = min(u(j-1, i), u(1, i));
            else
                b = min(u(j-1, i), u(j+1, i));
            end
            u(j,i) = local_update(a, b, f(j,i), h, u(j,i));
        end
    end
    % 4) sweep: i=Nx->1, j=Ny->1
    for jj = Ny:-1:1
        for ii = Nx:-1:1
            i = ii; j = jj;
            if boundary_mask(j,i); continue; end
            if(i==1)
                a = min(u(j, Nx), u(j, i+1));
            elseif(i==Nx)
                a = min(u(j, i-1), u(j, 1));
            else
                a = min(u(j, i-1), u(j, i+1));
            end

            if(j==1)
                b = min(u(Ny, i), u(j+1, i));

            elseif(j==Ny)
                b = min(u(j-1, i), u(1, i));
            else
                b = min(u(j-1, i), u(j+1, i));
            end
            u(j,i) = local_update(a, b, f(j,i), h, u(j,i));
        end
    end

    % error evaluation
    maxdiff = max(abs(u(:) - u_old(:)));
    fprintf('sweep %3d, maxdiff = %.3e\n', sweep, maxdiff);
    if maxdiff < tol
        fprintf('Converged after %d sweeps (tol=%.1e)\n', sweep, tol);
        break;
    end
end

if sweep == max_sweeps
    fprintf('Reached max_sweeps = %d (maxdiff = %.3e)\n', max_sweeps, maxdiff);
end

% for j=1:Ny
%     for i=1:Nx
%         if(phi(j,i)<0)
%             u(j,i)=-u(j,i);
%         end
%     end
% end

matrixphi=u;

%% visulization
figure('units','normalized','outerposition',[0 0 0.5 1]);
imagesc(x,y,u'); axis xy equal tight;
colorbar; hold on;
contour(x,y,u', 20, 'LineColor','k'); 
title('Initialization with Fast Sweeping Method');
plot(x(boundary_mask(1,:)), y(1, boundary_mask(1,:)), 'w.','MarkerSize',1); % dummy (not needed)
xlabel('x'); ylabel('y');

% overlay the exact boundary circle
th = linspace(0,2*pi,400);
plot(xc1 + r1*cos(th), yc1 + r1*sin(th), 'w--', 'LineWidth', 1.2);
plot(xc2 + r2*cos(th), yc2 + r2*sin(th), 'w--', 'LineWidth', 1.2);


%% end of main script

%% local_update from Zhao's paper
function u_new = local_update(a, b, f_ij, h, u_current)
    % ensure a <= b
    if a > b
        tmp = a; a = b; b = tmp;
    end
    % Delta = 2*(f*h)^2 - (a-b)^2
    Delta = 2*(f_ij * h)^2 - (a - b)^2;
    if Delta >= 0
        cand = (a + b + sqrt(max(0, Delta))) / 2;
    else
        cand = a + f_ij * h;
    end
    % upwind causality
    u_new = min(u_current, cand);
end

