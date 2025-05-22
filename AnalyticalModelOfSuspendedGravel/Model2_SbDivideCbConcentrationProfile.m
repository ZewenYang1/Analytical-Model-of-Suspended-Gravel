clc;
clear all;
close all;

% ==== 定义常量 ====
rho = 1000;     % 水密度 (kg/m3)
rhos = 2650;    % 泥沙密度 (kg/m3)
R = (rhos - rho) / rhos;  % 相对密度
A = 1.3e-7;     % Parker 参考浓度方程参数
g = 9.81;       % 重力加速度 (m/s?)
v = 1e-6;       % 动力粘度 (m?/s)
karmanconstant = 0.41; % 卡门常数
 beta =  1.0; % Schmidt 数
% beta =  0.49; % Schmidt 数

% ==== 水深与流速 ====
h = [50.42, 99.55, 93.12, 76.25, 54.63, 34.26, 18.80, 9.03];
U = [7.15, 10.05, 9.72, 8.79, 7.44, 5.89, 4.37, 3.02];

% ==== 粒径范围 ====
D = [0.001, 0.005, 0.010, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05];

% ==== 配色方案 ====
colors = flipud(parula(length(D)));

% ==== 创建图形窗口 ====
figure;
set(gcf, 'Position', [100, 100, 1400, 1000]);
rows = 2;
cols = 4;

% ==== 字体设置 ====
font_size = 20;
set(groot, 'DefaultAxesFontSize', font_size);
set(groot, 'DefaultTextFontSize', font_size);
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesXColor', 'k');
set(groot, 'DefaultAxesYColor', 'k');
set(groot, 'DefaultAxesZColor', 'k');
set(groot, 'DefaultTextColor', 'k');

% ==== 主循环 ====
for i = 1:length(h)
    u_shear = U(i) / 8.3;
    a = h(i) * 0.01;
    y = linspace(0, h(i), 100000);
    y_h = y / h(i);

    subplot(rows, cols, i);
    hold on;

    xtickformat('%.1f');
    ytickformat('%.1f');

    for j = 1:length(D)
        d = D(j) * ((R * g)^(1/3)) / (v^(2/3));
        Ws = 0.51 * (v / D(j)) * ((D(j)^3) * g * R / (v^2))^0.553;
        k = Ws / u_shear;
        P = k / (karmanconstant * beta);
        Cb = calculateCb(A, d, u_shear, Ws);  % 传入 u_shear 和 Ws
        Sb = calculateSb(y, h(i), a, Cb, P);

        % 绘图
        plot(Sb / Cb, y_h, 'Color', colors(j, :), 'LineWidth', 2.5);        
    end

    % 设置轴标签与图形属性
    xlabel('{\itS_b/C_b} (-)', 'FontSize', font_size, 'FontName', 'Arial');
    ylabel('{\ity/h} (-)', 'FontSize', font_size, 'FontName', 'Arial');
    set(gca, 'YDir', 'normal');
    box on;
    xlim([0 1]);
    ylim([0 1]);
    set(gca, 'TickLength', [0.03, 0.05]);
    hold off;
end

% ==== 添加 colorbar ====
c = colorbar('eastoutside');
caxis([min(D) max(D)]);
colormap(flipud(parula));
c.Label.String = 'Grain size {\itD} (mm)';
c.Label.FontSize = font_size;
c.Label.FontName = 'Arial';
tick_positions = linspace(min(D), max(D), length(D));
c.Ticks = tick_positions;
c.TickLabels = arrayfun(@(x) num2str(round(x * 1000)), D, 'UniformOutput', false);
c.Position = [0.92, 0.11, 0.018, 0.815];

% export_fig Beta=0.49_ConcentrationProfileDistribution.jpg -r1000 -transparent;

% ==== 函数定义 ====
function Cb = calculateCb(A, d, u_shear, Ws)
    Rf = 0.0738 * (d^1.021);
    ratio = u_shear / Ws;  
    Zu = 9.95 * ratio * (Rf^0.882);
    Cb = (A * Zu^5) / (1 + (A/0.3) * Zu^5);
end

function Sb = calculateSb(y, h1, a, Cb, P)
    Sb = (((h1 - y) ./ y) .* (a / (h1 - a)) .^ P) * Cb;
end


 
 
 

