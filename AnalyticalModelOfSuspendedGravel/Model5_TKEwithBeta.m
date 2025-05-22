clc;
clear all;
close all;

% ==== 基本常量 ====
g = 9.81;      % 重力加速度 (m/s^2)
rho = 1000;    % 水密度 (kg/m^3)
rhos = 2650;   % 颗粒密度 (kg/m^3)
R = (rhos - rho) / rho; 
v = 1e-6;      % 水的运动粘度 (m^2/s)
Mu = 0.001;
C_Nu = 0.09;

% ==== 卡门常数数组 ====
K_values = [0.41, 0.30];  % 两个不同的 K 值

% ==== 粗糙度高度 ====
z0 = 0.05 * 2.5;  % 固定粗糙度高度

% ==== 8 组水深与平均流速 ====
h_values = [50.42, 99.55, 93.12, 76.25, 54.63, 34.26, 18.80, 9.03];
U_values = [7.15, 10.05, 9.72, 8.79, 7.44, 5.89, 4.37, 3.02];

% ==== 图形窗口 ====
figure;
set(gcf, 'Position', [100, 90, 800, 1000]);

% 颜色数组用于区分不同 K
color_set = [255 145 37; 0 189 255] / 255;

% ==== 循环绘制子图 ====
for i = 1:8
    h = h_values(i);
    U = U_values(i);
    
    subplot(4, 2, i);
    hold on;

    % 存储绘图句柄
    plot_handles = gobjects(length(K_values), 1);

    for j = 1:length(K_values)
        K = K_values(j);

        % ==== 计算剪切流速 ====
        u_shear = U / 8.3;

        % ==== 计算 TKE ====
        y = linspace(0.05, h, 100);
        y_h = y / h;
        Epsilon = ((9.8 * exp(-3 * y / h)) ./ ((y / h) .^ 0.5)) .* ((u_shear ^ 3) / h);
        y1 = y * u_shear / v;
        u_rootmean = 2.3 * u_shear .* exp(-y / h) .* (1 - exp(-y1 / 10)) + 0.3 * y1 .* exp(-y1 / 10);
        Lambda = sqrt(15 * v * (u_rootmean.^2) ./ Epsilon);
        TKE = Epsilon .* K * h ./ ((C_Nu^0.5)*u_shear);

        % ==== 绘图并记录句柄 ====
        plot_handles(j) = plot(TKE, y_h, 'LineWidth', 2, 'Color', color_set(j,:));
    end

    % ==== 图形设置 ====
    set(gca, 'FontSize', 20);  
    xlabel('{\itTKE} (m^2/s^2)', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    ylabel('{\ity/h}', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    ytickformat('%.1f');
    xlim([0, 100]);
    set(gca, 'XTick', 0:20:100); 
    set(gca, 'YDir', 'normal');
    set(gca, 'TickLength', [0.03, 0.06]); 
    box on;

    % ==== 仅第一个子图添加图例 ====
    if i == 1
        legend(plot_handles, {'$\beta$ = 1', '$\beta$ = 0.73'}, 'Location', 'northeast', 'FontSize', 16, 'Interpreter', 'latex');
        %  K =0.41 is corresponds to β = 1; K =0.3 is corresponds to β = 0.73
        % (K/β always equal to 0.41, and β is a simple means to ‘keep K = 0.41 in the various equations. K is Karman constant)
                
        legend boxoff;   
    end
end

% export_fig TKE_with_β.jpg -r1000 -transparent;


