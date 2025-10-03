clc;
clear;
close all;

%% loading data
load('matrixphi_fourbubble.mat');
load('matrixF_fourbubble.mat');

%% main loop
[M,N]=size(matrixphi);
%count record the bubbles' tag
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


%% visulization
figure('Name','Bubble labels (Track) - discrete colors','NumberTitle','off','Renderer','painters');

maxLabel = max(Track(:));

bgColor = [1 1 1];      

if maxLabel == 0
    rgb = repmat(reshape(bgColor,1,1,3), M, N);
    image(rgb);
    axis equal tight off;
    title('No bubbles found (all zeros)');
    return;
end

nColorsNeeded = maxLabel;
baseColors = lines(max(7, nColorsNeeded)); 

if size(baseColors,1) < nColorsNeeded
    reps = ceil(nColorsNeeded / size(baseColors,1));
    baseColors = repmat(baseColors, reps, 1);
end
labelColors = baseColors(1:nColorsNeeded, :); 

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

bgmask = (Track == 0);
for c = 1:3
    ch = rgb(:,:,c);
    ch(bgmask) = bgColor(c);
    rgb(:,:,c) = ch;
end


image(rgb);
axis equal tight;
set(gca,'YDir','normal');  
title('CCL Bubble Tracking');

hold on;


props = regionprops(Track, 'Centroid', 'Area');
for k = 1:maxLabel
    mask = (Track == k);
    if ~any(mask(:)), continue; end

    B = bwboundaries(mask, 'noholes');
    for b = 1:numel(B)
        boundary = B{b};
        plot(boundary(:,2), boundary(:,1), 'k-', 'LineWidth', 0.7);
    end

    if k <= numel(props) && ~isempty(props(k).Centroid)
        c = props(k).Centroid;

        lblCol = labelColors(k,:);
        brightness = 0.299*lblCol(1) + 0.587*lblCol(2) + 0.114*lblCol(3);
        textColor = 'k';
        if brightness < 0.5, textColor = 'w'; end

        area_k = props(k).Area;
        fontSz = min(max(8, round(8 + sqrt(area_k)/6)), 18);

        text(c(1), c(2), num2str(k), 'Color', textColor, ...
            'FontWeight','bold', 'FontSize', fontSz, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
end

hold off;




%% recursion research

function S_tagged=search(matrixphi,j,i,flag)
persistent S;
persistent M;
persistent N;

if(flag==1)
    [M,N]=size(matrixphi);
    S=zeros(M,N);
end

S_tagged=[];

if(matrixphi(j,i)>0)
    return;
end

if(S(j,i)==1)
    return;
end


S(j,i)=1;


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

%Use flag to judge if it's the first recursion
if(flag==1)
    S_tagged=S;
end

end

%% add boundary cell

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

%% tagging function
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
