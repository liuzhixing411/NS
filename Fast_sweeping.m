clc;
clear;

% fast_sweep_demo.m
% Fast Sweeping demo using local analytic update (paper Eq. (2.4))
% Save as fast_sweep_demo.m and run in MATLAB.

clear; close all; clc;

%% 网格与参数
Lx = 8; Ly = 8;   % 物理域尺寸
Nx = 160; Ny = 160;   % 网格点数（越密越精确）
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
h = x(2)-x(1);
[X, Y] = meshgrid(x, y);

% 边界示例：圆形边界（半径 r0，圆周上设 u=0）
xc1 = 4; yc1 = 4;    % 圆心
r1 = 1;              % 半径

xc2=6;yc2=6;
r2=0.1;

% f(x)（slowness），这里取常数 1（即距离场）
f = ones(Ny, Nx);    % 注意 MATLAB 的矩阵索引 (row y, col x)

%% 初始化 u（未知点用 Inf）

u = Inf(Ny, Nx);

% 找到最接近圆周的格点，作为 Dirichlet 边界 u=0
phi = min(sqrt((X - xc1).^2 + (Y - yc1).^2)-r1,sqrt((X - xc2).^2 + (Y - yc2).^2)-r2);
boundary_mask = abs(phi) <= (h*0.6);  % 小带宽近似圆周
for j=1:Ny
    for i=1:Nx
        if(boundary_mask(j,i)==1)
            u(j,i)=phi(j,i);
        end
    end
end


%% Fast Sweeping 主循环（4 个方向交替）
tol = 1e-6;
max_sweeps = 100;   % 最多 sweep 次数
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

    % 收敛检测
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

%% 可视化 1：平面伪彩 + 等值线
figure('units','normalized','outerposition',[0 0 0.5 1]);
imagesc(x,y,u'); axis xy equal tight;
colorbar; hold on;
contour(x,y,u', 20, 'LineColor','k');  % 20 等值线
title('Initialization with Fast Sweeping Method');
plot(x(boundary_mask(1,:)), y(1, boundary_mask(1,:)), 'w.','MarkerSize',1); % dummy (not needed)
xlabel('x'); ylabel('y');

% overlay the exact boundary circle
th = linspace(0,2*pi,400);
plot(xc1 + r1*cos(th), yc1 + r1*sin(th), 'w--', 'LineWidth', 1.2);
plot(xc2 + r2*cos(th), yc2 + r2*sin(th), 'w--', 'LineWidth', 1.2);


%% end of main script

%% local_update 函数（实现论文式 (2.4) 的局部解析更新）
function u_new = local_update(a, b, f_ij, h, u_current)
    % 保证 a <= b
    if a > b
        tmp = a; a = b; b = tmp;
    end
    % 判别式 Delta = 2*(f*h)^2 - (a-b)^2
    Delta = 2*(f_ij * h)^2 - (a - b)^2;
    if Delta >= 0
        cand = (a + b + sqrt(max(0, Delta))) / 2;
    else
        cand = a + f_ij * h;
    end
    % 单调性保证（upwind causality）
    u_new = min(u_current, cand);
end
