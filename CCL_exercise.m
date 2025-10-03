clc;
clear;
close all;

%% 数据导入
load('matrixphi_fourbubble.mat');
load('matrixF_fourbubble.mat');

%% 主函数
[M,N]=size(matrixphi);
%count记录当前气泡标号
count=0;
Track=zeros(M,N);

for j=1:M
    for i=1:N
        if(matrixphi(j,i)<0 && Track(j,i)==0)
            count=count+1;
            S_tagged=search(matrixphi,j,i,1);
            Gridlabel=Add(S_tagged,matrixF,matrixphi);
            Track=tag(Gridlabel,count,Track);
        end
    end
end


%% 可视化（离散颜色 + 背景单独颜色 + 保留编号）
figure('Name','Bubble labels (Track) - discrete colors','NumberTitle','off','Renderer','painters');

% 基本信息
maxLabel = max(Track(:));

% 设置背景颜色（RGB），例如浅灰或白色
bgColor = [1 1 1];      % 白色背景。改为 [0.95 0.95 0.95] 等也行

% 若无气泡，直接显示背景并退出
if maxLabel == 0
    rgb = repmat(reshape(bgColor,1,1,3), M, N);
    image(rgb);
    axis equal tight off;
    title('No bubbles found (all zeros)');
    return;
end

% 为每个标签生成离散颜色（不包含背景）
% 使用 lines / hsv / parula 等生成 distinct colors
nColorsNeeded = maxLabel;
baseColors = lines(max(7, nColorsNeeded));  % lines 常给较区分颜色；确保最少7色以便更好区分

% 如果生成的颜色少于需要，则循环重复（通常不会发生）
if size(baseColors,1) < nColorsNeeded
    reps = ceil(nColorsNeeded / size(baseColors,1));
    baseColors = repmat(baseColors, reps, 1);
end
labelColors = baseColors(1:nColorsNeeded, :);  % 每个标签对应一行 RGB

% 构造 RGB 图像（M x N x 3）
rgb = zeros(M, N, 3);
for k = 1:maxLabel
    mask = (Track == k);
    if any(mask(:))
        for c = 1:3
            channel = rgb(:,:,c);
            channel(mask) = labelColors(k,c);
            rgb(:,:,c) = channel;
        end
    end
end
% 背景色填充 (Track == 0)
bgmask = (Track == 0);
for c = 1:3
    ch = rgb(:,:,c);
    ch(bgmask) = bgColor(c);
    rgb(:,:,c) = ch;
end

% 显示图像
image(rgb);
axis equal tight;
set(gca,'YDir','normal');  % 保持行列方向直观
title('CCL Bubble Tracking');

hold on;

% 画轮廓 & 标注编号（质心）
props = regionprops(Track, 'Centroid', 'Area');
for k = 1:maxLabel
    mask = (Track == k);
    if ~any(mask(:)), continue; end

    % 轮廓（多连通片时多条边界）
    B = bwboundaries(mask, 'noholes');
    for b = 1:numel(B)
        boundary = B{b};
        % 画边界：用黑色或白色边线以保证对比
        plot(boundary(:,2), boundary(:,1), 'k-', 'LineWidth', 0.7);
    end

    % 标注质心（regionprops 返回 Centroid = [x y] = [col row]）
    if k <= numel(props) && ~isempty(props(k).Centroid)
        c = props(k).Centroid;
        % 根据该标签颜色亮度选择文本颜色（黑/白）
        lblCol = labelColors(k,:);
        brightness = 0.299*lblCol(1) + 0.587*lblCol(2) + 0.114*lblCol(3);
        textColor = 'k';
        if brightness < 0.5, textColor = 'w'; end

        % 字体大小与区域大小略相关，避免对大/小区域字体过大或过小
        area_k = props(k).Area;
        fontSz = min(max(8, round(8 + sqrt(area_k)/6)), 18);

        text(c(1), c(2), num2str(k), 'Color', textColor, ...
            'FontWeight','bold', 'FontSize', fontSz, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end

hold off;




%% 写一个递归搜索的函数

function S_tagged=search(matrixphi,j,i,flag)
%S对访问过的位置标记1，其余为0
persistent S;
persistent M;
persistent N;

if(flag==1)
    [M,N]=size(matrixphi);
    S=zeros(M,N);
end

S_tagged=[];



%若距离为正也跳出
if(matrixphi(j,i)>0)
    return;
end

%若标记过则跳出
if(S(j,i)==1)
    return;
end

%标记
S(j,i)=1;

%递归,上下左右
if(j==M)
    [~]=search(matrixphi,1,i,2);
else
    [~]=search(matrixphi,j+1,i,2);
end

if(j==1)
    [~]=search(matrixphi,M,i,2);
else
    [~]=search(matrixphi,j-1,i,2);
end

if(i==1)
    [~]=search(matrixphi,j,N,2);

else
    [~]=search(matrixphi,j,i-1,2);
end

if(i==N)
    [~]=search(matrixphi,j,1,2);

else
    [~]=search(matrixphi,j,i+1,2);
end

%用flag判断是不是第一层
if(flag==1)
    S_tagged=S;
end

end

%% 写一个找补充边界点的函数

function Gridlabel=Add(S,matrixF,matrixphi)
[M,N]=size(S);
Gridlabel=S;

for j=1:M
    for i=1:N
        if(S(j,i)==1)

            if(j==M)
                if(matrixF(1,i)<1 && matrixphi(1,i)>0)
                    Gridlabel(1,i)=1;
                end
            else

                if(matrixF(j+1,i)<1 && matrixphi(j+1,i)>0)
                    Gridlabel(j+1,i)=1;
                end

            end

            if(j==1)
                if(matrixF(M,i)<1 && matrixphi(M,i)>0)
                    Gridlabel(M,i)=1;
                end
            else

                if(matrixF(j-1,i)<1 && matrixphi(j-1,i)>0)
                    Gridlabel(j-1,i)=1;
                end

            end

            if(i==N)
                if(matrixF(j,1)<1 && matrixphi(j,1)>0)
                    Gridlabel(j,1)=1;
                end
            else

                if(matrixF(j,i+1)<1 && matrixphi(j,i+1)>0)
                    Gridlabel(j,i+1)=1;
                end

            end

            if(i==1)
                if(matrixF(j,N)<1 && matrixphi(j,N)>0)
                    Gridlabel(j,N)=1;
                end
            else

                if(matrixF(j,i-1)<1 && matrixphi(j,i-1)>0)
                    Gridlabel(j,i-1)=1;
                end

            end

        end
    end
end

end

%% 写一个标记函数
function Track=tag(Gridlabel,tagnumber,track)

[M,N]=size(track);
Track=track;
for j=1:M
    for i=1:N
        if(Gridlabel(j,i)==1)
            Track(j,i)=tagnumber;
        end
    end
end
end