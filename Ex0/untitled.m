clc;
clear;
y_start = 0;
y_end = 1;
formatSpec = '%f';
fileID = fopen("output_D_double_10.txt", "r");
y = fscanf(fileID,formatSpec);
fclose(fileID);

for i = 1:length(y)
    x(i) = 0+i*(y_end-y_start)/(length(y)+1);
end

plot(x, y,'LineWidth',2,'Color',"#0072BD")
plot(x,x+y(1)-x(1))
plot(x,-1.*x+1*x(end)+y(end))
%%
formatSpec = '%f';
y_start = 0;
y_end = 1;

fig1 = figure ("Name","Solution For The Given Equation",'Position',[100 100 900 500]);
subplot(1,2,1)
    axis equal
    hold all
    
    fileID = fopen("output_D_double_5.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,'-.','LineWidth',2,'Color',"#0072BD")
% ------------------------------------------------------
    fileID = fopen("output_D_double_10.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,'LineWidth',1.75,'Color',"#D95319")
    % ------------------------------------------------------
    fileID = fopen("output_D_double_100.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,'--','LineWidth',1.5,'Color',"#EDB120")
% ------------------------------------------------------
    fileID = fopen("output_D_double_1000.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,':','LineWidth',1.25,'Color',"#77AC30")

    grid on
    grid minor
    ylabel("Y [-]")
    xlabel("X [-]")
    title("Solution For Dirichlet Boundary Conditions For Different N; Y_0 = 0, Y_N = 1; double")
    subtitle("Almog Dobrescu 214254252")
    legend({'N = 5','N = 10','N = 100','N = 1000'},'FontSize',11 ,'Location','southeast')
% ########################################################################
subplot(1,2,2)
    axis equal
    hold all
    
    fileID = fopen("output_D_float_5.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,'-.','LineWidth',2,'Color',"#0072BD")
% ------------------------------------------------------
    fileID = fopen("output_D_float_10.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,'LineWidth',1.75,'Color',"#D95319")
    % ------------------------------------------------------
    fileID = fopen("output_D_float_100.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,'--','LineWidth',1.5,'Color',"#EDB120")
% ------------------------------------------------------
    fileID = fopen("output_D_float_1000.txt", "r");
    y = [];
    x = [];
    y1 = [];
    x1 = [];
    y = fscanf(fileID,formatSpec);
    for i = 1:length(y)
        x(i) = 0+i*(1-0)/(length(y)+1);
    end
    fclose(fileID);
    
    for i = 1:length(y)
        y1(i+1) = y(i);
        x1(i+1) = x(i);
    end
    
    y1(1) = y_start;
    x1(1) = 0;
    y1(end+1) = y_end;
    x1(end+1) = 1;
    
    plot(x1, y1,':','LineWidth',1.25,'Color',"#77AC30")

    grid on
    grid minor
    ylabel("Y [-]")
    xlabel("X [-]")
    title("Solution For Dirichlet Boundary Conditions For Different N; Y_0 = 0, Y_N = 1; float")
    subtitle("Almog Dobrescu 214254252")
    legend({'N = 5','N = 10','N = 100','N = 1000'},'FontSize',11 ,'Location','southeast')
%exportgraphics(fig1, 'grap1.png','Resolution',1800);

%% 
formatSpec = '%f';
y_start = 0;
y_end = 1;

fig2 = figure ("Name","Solution For The Given Equation",'Position',[300 100 900 500]);
subplot(1,2,1)
    axis equal
    hold all
    
    fileID = fopen("output_N_double_5.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#0072BD")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#0072BD", 'HandleVisibility', 'off')    
    plot(x, y,'-.','LineWidth',2,'Color',"#0072BD")
% ------------------------------------------------------
    fileID = fopen("output_N_double_10.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#D95319")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#D95319", 'HandleVisibility', 'off')    
    plot(x, y,'LineWidth',1.75,'Color',"#D95319")
    % ------------------------------------------------------
    fileID = fopen("output_N_double_100.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#EDB120")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#EDB120", 'HandleVisibility', 'off')    
    plot(x, y,'-','LineWidth',1.5,'Color',"#EDB120")
% ------------------------------------------------------
    fileID = fopen("output_N_double_1000.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#77AC30")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#77AC30", 'HandleVisibility', 'off')    
    plot(x, y,'-.','LineWidth',1.5,'Color',"#77AC30")

    grid on
    grid minor
    ylabel("Y [-]")
    xlabel("X [-]")
    title("Solution For Neumann Boundary Conditions For Different N; Y'_0 = 1, Y'_N = -1; double")
    subtitle("Almog Dobrescu 214254252")
    legend({'y=+/-x+b','N = 5','y=+/-x+b','N = 10','y=+/-x+b','N = 100','y=+/-x+b','N = 1000'},'FontSize',11 ,'Location','southeast')
    ylim([4.9 6.4]) 
% ########################################################################
subplot(1,2,2)
axis equal
    hold all
    
    fileID = fopen("output_N_float_5.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#0072BD")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#0072BD", 'HandleVisibility', 'off')    
    plot(x, y,'-.','LineWidth',2,'Color',"#0072BD")
% ------------------------------------------------------
    fileID = fopen("output_N_float_10.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#D95319")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#D95319", 'HandleVisibility', 'off')    
    plot(x, y,'LineWidth',1.75,'Color',"#D95319")
    % ------------------------------------------------------
    fileID = fopen("output_N_float_100.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#EDB120")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#EDB120", 'HandleVisibility', 'off')    
    plot(x, y,'-','LineWidth',1.5,'Color',"#EDB120")
% ------------------------------------------------------
    fileID = fopen("output_N_float_1000.txt", "r");
    y = [];
    x = [];
    y = fscanf(fileID,formatSpec);
    x = linspace(0, 1, length(y));
    fclose(fileID);
    
    plot(x,x+y(1)-x(1),'-.','LineWidth',1,'Color',"#77AC30")
    plot(x,-1.*x+1*x(end)+y(end),'-.','LineWidth',1,'Color',"#77AC30", 'HandleVisibility', 'off')    
    plot(x, y,'-.','LineWidth',1.5,'Color',"#77AC30")

    grid on
    grid minor
    ylabel("Y [-]")
    xlabel("X [-]")
    title("Solution For Neumann Boundary Conditions For Different N; Y'_0 = 1, Y'_N = -1; float")
    subtitle("Almog Dobrescu 214254252")
    legend({'y=+/-x+b','N = 5','y=+/-x+b','N = 10','y=+/-x+b','N = 100','y=+/-x+b','N = 1000'},'FontSize',11 ,'Location','southeast')
    ylim([4.9 6.4]) 
%exportgraphics(fig2, 'grap2.png','Resolution',1000);

%%
formatSpec = '%f';

fig3 = figure ("Name","Solution For The Given Equation",'Position',[500 100 900 500]);
subplot(2,1,1)
fileID = fopen("output_D_float_1000_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#0072BD")
hold on
fileID = fopen("output_D_float_100_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#EDB120")
hold on
fileID = fopen("output_D_double_1000_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#D95319")
hold on
fileID = fopen("output_D_double_100_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#77AC30")


grid on
grid minor
ylabel("Y [-]")
xlabel("X [-]")
title("Dirichlet - Absolut Value of The Error For Double And Float For Different N")
subtitle("Almog Dobrescu 214254252")
legend({'float, N = 1000','float, N = 100','double, N = 1000','double, N = 100'},'FontSize',11 ,'Location','southeast')

subplot(2,1,2)
fileID = fopen("output_N_float_1000_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#0072BD")
hold on
fileID = fopen("output_N_float_100_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#EDB120")
hold on
fileID = fopen("output_N_double_1000_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#D95319")
hold on
fileID = fopen("output_N_double_100_check.txt", "r");
y = [];
x = [];
y = fscanf(fileID,formatSpec);
x = linspace(0, 1, length(y));
fclose(fileID);
semilogy(x, abs(y),'','LineWidth',1.5,'Color',"#77AC30")

grid on
grid minor
ylabel("Y [-]")
xlabel("X [-]")
title("Neumann - Absolut Value of The Error For Double And Float For Different N")
subtitle("Almog Dobrescu 214254252")
legend({'float, N = 1000','float, N = 100','double, N = 1000','double, N = 100'},'FontSize',11 ,'Location','southeast')
%exportgraphics(fig3, 'grap3.png','Resolution',1000);

