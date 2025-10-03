clc;
clear;

%% 数据导入
load('matrixF.mat');
load('matrixphi.mat');

%% 基本参数设置
%连续连通域不被删除的最小格子数
maxSize=50;
epsilon = 0.02;
h=0.05;
%% 主代码

%找到需要更改区域
groups_lin=small_zero_groups(matrixF,maxSize);

kink_mask = find_kinks(matrixphi, epsilon, h);

% groups_lin 已由 small_zero_groups(matrixF,maxSize) 得到
% groups_lin 是 cell array

% 1) 找出所有封闭环（并筛选出包含 groups_lin 的那些）
regions = find_kink_regions(kink_mask, groups_lin);

% 2) 可视化（把 groups_lin 传进去用于高亮）
plot_kink_regions(matrixphi, kink_mask, regions, groups_lin);

%% 函数
function kink_mask = find_kinks(matrixphi, epsilon, h)
% FIND_KINKS  找出 phi 矩阵中的 kink 节点（周期性 BC）
%
% Usage:
%   kink_mask = find_kinks(phi, epsilon)
%   kink_mask = find_kinks(phi, epsilon, h)
%
% Inputs:
%   phi      - MxN 矩阵（你的 160x160 测试矩阵）
%   epsilon  - 阈值（比如 0.01）
%   h        - 网格间距，optional（默认 1）
%
% Output:
%   kink_mask - MxN logical 矩阵，true 表示该节点为 kink
%
% 例：
%   kink_mask = find_kinks(phi, 0.02);
%   regions = find_kink_regions(kink_mask);
%   plot_kink_regions(phi, kink_mask, regions);

    if nargin < 3 || isempty(h), h = 1; end
    if nargin < 2
        error('Need phi and epsilon.');
    end

    % 使用周期性差分（circshift），注意 MATLAB 矩阵 (row=y, col=x)
    % d/dx: shift left/right on columns, d/dy: shift up/down on rows
    dphi_dx = (circshift(matrixphi, [0,-1]) - circshift(matrixphi, [0,1])) / (2*h);
    dphi_dy = (circshift(matrixphi, [-1,0]) - circshift(matrixphi, [1,0])) / (2*h);

    gradmag = sqrt(dphi_dx.^2 + dphi_dy.^2);   % >=0
    % 判断 kink
    kink_mask = abs(1 - abs(gradmag)) > epsilon;
    kink_mask = logical(kink_mask);

end


%% 写一个找删除域的函数
function groups_lin = small_zero_groups(A, maxSize)
    % 只关注不等于1的位置
    bw = (A ~= 1);
    % 4-连通
    CC = bwconncomp(bw, 4);
    sizes = cellfun(@numel, CC.PixelIdxList);
    keep = find(sizes < maxSize);
    groups_lin = CC.PixelIdxList(keep);
end

function regions = find_kink_regions(kink_mask, groups_lin)
    if nargin < 2
        groups_lin = {};
    end

    [M, N] = size(kink_mask);
    [r_list, c_list] = find(kink_mask);
    num_nodes = numel(r_list);
    regions = {};

    if num_nodes == 0
        return;
    end

    % 预处理：筛选出不在kink_mask上的groups
    valid_groups = {};
    for g = 1:numel(groups_lin)
        linlist = groups_lin{g};
        if ~any(kink_mask(linlist)) % 确保group中没有点在kink_mask上
            valid_groups{end+1} = linlist;
        end
    end
    
    % 如果没有有效的groups，则返回所有环
    if isempty(valid_groups)
        return_all = true;
    else
        return_all = false;
    end

    % 创建线性索引到节点索引的映射
    lin_idx = sub2ind([M,N], r_list, c_list);
    lin2node = containers.Map('KeyType','uint32','ValueType','uint32');
    for k = 1:num_nodes
        lin2node(uint32(lin_idx(k))) = uint32(k);
    end

    % 构建邻接矩阵（省略部分代码保持不变）
    % ... [保持原有的邻接矩阵构建代码]
     % 构造无向稀疏邻接矩阵（对称）
    % 预分配邻接边 lists
    I = []; J = [];

    % 8 邻域偏移（row, col）
    nbrs = [ -1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1 ];

    for k = 1:num_nodes
        r = r_list(k); c = c_list(k);
        for t = 1:size(nbrs,1)
            rr = r + nbrs(t,1);
            cc = c + nbrs(t,2);
            % 周期 wrap
            if rr < 1, rr = rr + M; end
            if rr > M, rr = rr - M; end
            if cc < 1, cc = cc + N; end
            if cc > N, cc = cc - N; end
            % 若邻点也是 kink，则加入边
            if kink_mask(rr, cc)
                lin_n = uint32(sub2ind([M,N], rr, cc));
                if lin2node.isKey(lin_n)
                    j = lin2node(lin_n);
                    I(end+1) = k; %#ok<AGROW>
                    J(end+1) = j; %#ok<AGROW>
                end
            end
        end
    end

    % 构造无向稀疏邻接矩阵（对称）
    A = sparse([I J], [J I], 1, num_nodes, num_nodes);

    % 找连通分量
    G = graph(A);
    comps = conncomp(G);
    comp_ids = unique(comps);

    region_count = 0;
    for cid = comp_ids
        nodes_in_comp = find(comps == cid);
        ncomp = numel(nodes_in_comp);
        if ncomp <= 2
            continue;
        end
        
        % 检查是否为环（所有节点度数为2）
        subA = A(nodes_in_comp, nodes_in_comp);
        degs = sum(subA,2);
        if ~all(degs == 2)
            continue;
        end

        % 提取环的坐标（省略顺序化代码）
        % ... [保持原有的环提取代码]
        
        % 解包环坐标以适应周期性
        [RRu, CCu] = unwrap_polygon_coords(RR, CC, M, N);
        xpoly = CCu(:);
        ypoly = RRu(:);

        % 计算环的内部点
        medC = median(xpoly);
        medR = median(ypoly);
        
        [Cgrid, Rgrid] = meshgrid(1:N, 1:M);
        Cg = double(Cgrid(:)); 
        Rg = double(Rgrid(:));
        
        Cg_shift = shift_coords_to_center(Cg, medC, N);
        Rg_shift = shift_coords_to_center(Rg, medR, M);
        
        [in_all, on_all] = inpolygon(Cg_shift, Rg_shift, xpoly, ypoly);
        in_all = reshape(in_all, M, N);
        on_all = reshape(on_all, M, N);
        interior_mask = in_all & ~on_all;

        % 检查环是否包含所有有效groups
        accept = return_all; % 如果没有有效groups，接受所有环
        
        if ~return_all
            accept = true;
            for g = 1:numel(valid_groups)
                linlist = valid_groups{g};
                [rg, cg] = ind2sub([M,N], linlist);
                
                % 检查group中的所有点是否都在环内部
                for i = 1:numel(linlist)
                    if ~interior_mask(rg(i), cg(i))
                        accept = false;
                        break;
                    end
                end
                
                if ~accept
                    break;
                end
            end
        end

        if accept
            region_count = region_count + 1;
            boundary = [RR(:), CC(:)];
            [ir, ic] = find(interior_mask);
            interior = [ir, ic];
            allcoords = unique([boundary; interior], 'rows', 'stable');
            regions{region_count} = struct('boundary', boundary, ...
                                          'interior', interior, ...
                                          'all', allcoords);
        end
    end

    if isempty(regions)
        regions = {};
    end
end


%% ---------- Helper: unwrap polygon coords (recentre around median) ----------
function [Rshift, Cshift] = unwrap_polygon_coords(R, C, M, N)
    % 将行列坐标 R,C recentre 到一个连续区间（避免多段跨域）
    % 基本策略：以列中值 median(C) 为中心，把所有列移动到距中心不超过 N/2 的同一周期
    medC = median(C);
    medR = median(R);
    Cshift = shift_coords_to_center(C, medC, N);
    Rshift = shift_coords_to_center(R, medR, M);
end

%% ---------- Helper: shift_coords_to_center ----------
function coords_shift = shift_coords_to_center(coords, center, period)
    % coords: vector of coords (double)
    % center: scalar center (double)
    % period: period (N for columns or M for rows)
    coords_shift = coords;
    % while difference too large, shift by +/- period
    % use vectorized approach: compute delta = coords - center, then add/sub multiples
    delta = coords_shift - center;
    % compute integer number of periods to add so that |delta - k*period| <= period/2
    k = round(delta / period);
    coords_shift = coords_shift - k * period;
end

%% 替换后的 plot_kink_regions（绘制边界与内部点并高亮 groups）
function plot_kink_regions(matrixphi, kink_mask, regions, groups_lin)
% PLOT_KINK_REGIONS 可视化 phi、kink 点与被接受的封闭环（并展示环内部）
% Usage:
%   plot_kink_regions(phi, kink_mask, regions)
%   plot_kink_regions(phi, kink_mask, regions, groups_lin)
    if nargin < 4, groups_lin = {}; end
    figure; imagesc(matrixphi); axis xy equal tight; colormap parula; colorbar;
    title('phi with kink points and detected closed kink regions');
    hold on;
    [r_all, c_all] = find(kink_mask);
    plot(c_all, r_all, '.k', 'MarkerSize', 6); % kink 点（黑点）

    if isempty(regions)
        legend('kink points');
        hold off;
        return;
    end

    nreg = numel(regions);
    cmap = hsv(max(3,nreg));
    for k = 1:nreg
        s = regions{k};
        if isempty(s), continue; end
        % boundary
        b = s.boundary;
        rr = b(:,1); cc = b(:,2);
        rr_closed = [rr; rr(1)]; cc_closed = [cc; cc(1)];
        plot(cc_closed, rr_closed, '-', 'LineWidth', 2, 'Color', cmap(mod(k-1,size(cmap,1))+1,:));
        % interior points as small dots (same color but lighter)
        inpts = s.interior;
        if ~isempty(inpts)
            plot(inpts(:,2), inpts(:,1), '.', 'MarkerSize', 6, 'Color', cmap(mod(k-1,size(cmap,1))+1,:));
        end
        % mark start of boundary
        plot(cc(1), rr(1), 'o', 'MarkerFaceColor', cmap(mod(k-1,size(cmap,1))+1,:), 'MarkerEdgeColor','k');
    end

    % 如果提供 groups_lin，把这些 group 的质心标注出来以便检查
    for g = 1:numel(groups_lin)
        linlist = groups_lin{g};
        if isempty(linlist), continue; end
        [rg, cg] = ind2sub(size(matrixphi), linlist);
        % group centroid
        plot(mean(cg), mean(rg), 'ks', 'MarkerSize',10, 'MarkerFaceColor','y');
    end

    legend('kink points', 'closed region boundaries', 'region interior dots');
    hold off;
end






