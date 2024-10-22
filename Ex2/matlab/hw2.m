clc; clear; close all;
ni = 51;
nj = 26;
gamma = 1.4;
formatSpec = '%f';

fileID = fopen("results/Q0_mat.txt", "r");
Q0s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q0(i,:) = Q0s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q0;

fileID = fopen("results/Q1_mat.txt", "r");
Q1s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q1(i,:) = Q1s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q1;

fileID = fopen("results/Q2_mat.txt", "r");
Q2s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q2(i,:) = Q2s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q2;

fileID = fopen("results/Q3_mat.txt", "r");
Q3s = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    Q3(i,:) = Q3s(1+(i-1)*ni:ni+(i-1)*ni);
end
Q3;

fileID = fopen("results/x_mat.txt", "r");
xs = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
   x(i,:) = xs(1+(i-1)*ni:ni+(i-1)*ni);
end
x;

fileID = fopen("results/y_mat.txt", "r");
ys = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    y(i,:) = ys(1+(i-1)*ni:ni+(i-1)*ni);
end
y;

fileID = fopen("results/U_mat.txt", "r");
Us = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    U(i,:) = Us(1+(i-1)*ni:ni+(i-1)*ni);
end
U;

fileID = fopen("results/V_mat.txt", "r");
Vs = fscanf(fileID,formatSpec);
fclose(fileID);
for i = 1:nj
    V(i,:) = Vs(1+(i-1)*ni:ni+(i-1)*ni);
end
V;

% -------------------------------------------

u = Q1./Q0;
v = Q2./Q0;
p = (gamma - 1).*(Q3 - 0.5.*Q0.*(u.^2 + v.^2));
a = sqrt(gamma.*p./Q0);
M = sqrt(u.^2 + v.^2)./a;

fig1 = figure ("Name","1",'Position',[500 100 900 500]);
axis equal
hold all
contourf(x, y, M, 2000, "LineStyle","none");
colorbar;
colormap("jet");

scale = 3/10000;
q = quiver(x, y, u*scale, v*scale, "off", "Color", "#FFFFFF");
xlim([-0.25, 1.1])
ylim([-0.3, 0.3])

plot(x(:,end), y(:,end),'-*','Color',"#7E2F8E")
plot(x(:,1), y(:,1),'-*m')
plot(x(end,:), y(end,:),'-*r')
plot(x(1,:), y(1,:),'-k')
% 
% for i = 2:(length(x(1,:))-1)
%     plot(x(:,i), y(:,i),'-','Color',"#0072BD")
% end
% for i = 2:(length(x(:,1))-1)
%     plot(x(i,:), y(i,:),'-','Color',"#0072BD")
% end

%%

fig2 = figure ("Name","2",'Position',[100 350 900 500]);
% [xstart, ystart] = meshgrid(25:0.075:27, nj);
[xstart, ystart] = meshgrid(((ni+1)/2-1):0.075:((ni+1)/2+1), nj);
quiver(U, V)
sl1 = stream2(U, V, xstart, ystart,[0.0005, 4e5]);
line = streamline(sl1);
for i = 1:length(line)
    line(i).Color = "r";
end
close name 2

sl = sl1;

size_sl = size(sl);
for k = 1:size_sl(1)
    k;
    for j = 1:length(sl)
        x_val = [];
        y_val = [];
        for i = 1:length(sl{k,j}(:, 1))
            
            x_index_low = floor(sl{k,j}(i, 1));
            x_index = sl{k,j}(i, 1);
            x_index_high = ceil(sl{k,j}(i, 1));
        
            y_index_low = floor(sl{k,j}(i, 2));
            y_index = sl{k,j}(i, 2);
            y_index_high = ceil(sl{k,j}(i, 2));
            if (x_index_high == x_index_low)
                x_val_y_low = x(y_index_low, x_index);
                x_val_y_high = x(y_index_high, x_index);
        
                y_val_y_low = y(y_index_low, x_index);
                y_val_y_high = y(y_index_high, x_index);
            else
                x_val_y_low = x(y_index_low, x_index_low) + (x(y_index_low, x_index_low) - x(y_index_low, x_index_high)) / (x_index_low - x_index_high) * (x_index - x_index_low);
                x_val_y_high = x(y_index_high, x_index_low) + (x(y_index_high, x_index_low) - x(y_index_high, x_index_high)) / (x_index_low - x_index_high) * (x_index - x_index_low);
        
                y_val_y_low = y(y_index_low, x_index_low) + (y(y_index_low, x_index_low) - y(y_index_low, x_index_high)) / (x_index_low - x_index_high) * (x_index - x_index_low);
                y_val_y_high = y(y_index_high, x_index_low) + (y(y_index_high, x_index_low) - y(y_index_high, x_index_high)) / (x_index_low - x_index_high) * (x_index - x_index_low);
            end
            if (y_index_high == y_index_low)
                x_val(i) = x_val_y_low;
        
                y_val(i) = y_val_y_low;
            else
                x_val(i) = x_val_y_low + (x_val_y_low - x_val_y_high) / (y_index_low - y_index_high) * (y_index - y_index_low);
                
                y_val(i) = y_val_y_low + (y_val_y_low - y_val_y_high) / (y_index_low - y_index_high) * (y_index - y_index_low);
            end
        end
        
        figure(1)
        hold all
        if k == 1
            plot (x_val, y_val, "--", "LineWidth", 1, "Color", "k")
        end
    end
end