clc;
clear all;
close all;

% ==== 常量定义 ====
rho = 1000;     % 水密度 (kg/m^3)
rhos = 2650;    % 泥沙密度 (kg/m^3)
R = (rhos - rho) / rho;  % 相对密度
A = 1.3e-7;     % Parker参考浓度方程参数
g = 9.81;       % 重力加速度 (m/s^2)
v = 1e-6;       % 水的运动黏度 (m^2/s)
karmanconstant = 0.41;   % 卡门常数

% ==== 输入数据 ====    
h = [50.42, 99.55, 93.12, 76.25, 54.63, 34.26, 18.80, 9.03];% 水深 (m)
U = [7.15, 10.05, 9.72, 8.79, 7.44, 5.89, 4.37, 3.02];% 水深 (m)
D = [0.001, 0.005, 0.010, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]; % 粒径 (m)

% ==== 初始化结果矩阵 ====
C_avg_matrix = zeros(length(h), length(D));
ratio1_matrix = zeros(length(h), length(D));

% ==== 计算 C_avg 和 G_star/R_star^2 ====
for i = 1:length(h)
    u_shear = U(i) / 8.3;
    a = h(i) * 0.01;
    y = linspace(0, h(i), 10000);

    for j = 1:length(D)
        D_now = D(j);
        d = D_now * ((R * g)^(1/3)) / (v^(2/3));
        Ws = 0.51 * (v / D_now) * ((D_now^3 * g * R) / (v^2))^0.553;
        k_mov = Ws / u_shear;
        P = k_mov / karmanconstant;

        % 参考浓度
        Cb = calculateCb(A, d);

        % Sb 剖面与积分
        y = linspace(a, h(i), 10000); % 从 参考高度a积分到水面
        Sb = calculateSb(y, h(i), a, Cb, P);
        Sb_fun = @(yy) interp1(y, Sb, yy, 'pchip');
        C_avg = integral(Sb_fun, a, h(i), 'RelTol', 1e-8, 'AbsTol', 1e-12) / (h(i) - a);
        C_avg_matrix(i, j) = C_avg*100;

        % G_star / R_star^2 计算
        z = a * u_shear / v;
        u_rms = u_shear * 2.3 * exp(-a/h(i)) * (1 - exp(-z / 10)) + 0.3 * z * exp(-z / 10);
        Epsilon = (9.8 * exp(-3 * a / h(i)) * u_shear^3) / ((a / h(i))^0.5 * h(i));
        Lambda = sqrt(15 * v * u_rms^2 / Epsilon);
        R_star = Ws / u_rms;
        G_star = g * R * Lambda * Cb / (u_rms^2);
        ratio1_matrix(i, j) = G_star / R_star^2;
    end
end

% ==== 图：8个子图展示 C_avg vs G_star/R_star^2 ====
figure;
set(gcf, 'Position', [200, 100, 1600, 800]);
rows = 2; cols = 4;

% ==== 颜色映射 ====
nD = length(D);  % 粒径种类数
map = mymap('viridis');  % 自定义 colormap
cmap = map(round(linspace(1, size(map,1), nD)), :);  % 为每个粒径分配颜色
tick_vals = D * 1000;

for i = 1:length(h)
    subplot(rows, cols, i);
    xdata = ratio1_matrix(i, :);
    ydata = C_avg_matrix(i, :);
    valid = isfinite(xdata) & isfinite(ydata);
    scatter(xdata(valid), ydata(valid), 100, tick_vals(valid), 'filled', 'MarkerEdgeColor', 'k');
    hold on;
    if i <= 5
        ylim([1e-8, 1e4]);
        set(gca, 'YTick', [1e-8 1e-4 1e0 1e4]);
        xlim([1e-3, 1e1]);
        set(gca, 'XTick', [1e-3 1e-1 1e1]);
    else
        ylim([1e-20, 1e2]);
        set(gca, 'YTick', [1e-20 1e-13 1e-6 1e2]);
        xlim([1e-3, 1e1]);
        set(gca, 'XTick', [1e-3 1e-1 1e1]);
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('{\itG/R^{2}}','FontSize',20);
    ylabel('{\itC^*} (%)','FontSize',20);
    set(gca, 'TickLength', [0.04, 0.06]);  
    colormap(map);
    caxis([min(tick_vals), max(tick_vals)]);
    xRange = get(gca, 'XLim');
    plot(xRange, [10 10], '-', 'Color',[0.7 0.7 0.7],'LineWidth', 2.5);
    plot(xRange, [1 1],  '-', 'Color',[0.8 0.9 1.0], 'LineWidth', 2.5);
    scatter(xdata(valid), ydata(valid), 100, tick_vals(valid), 'filled', 'MarkerEdgeColor', 'k');
    scatter(xdata(valid), ydata(valid), 100, tick_vals(valid), 'filled', 'MarkerEdgeColor', 'k');
    box on;
end

% ==== 添加 colorbar ====
hcb = colorbar('Position', [0.93 0.141 0.0187 0.782]);
hcb.Label.String = 'Grain size {\itD} (mm)';
hcb.Label.FontSize = 20;
hcb.FontSize = 20;
hcb.Ticks = tick_vals;
hcb.TickLabels = arrayfun(@(x) num2str(x), tick_vals, 'UniformOutput', false);

export_fig C_starWithGR2.jpg -r1000 -transparent;

%==== 函数定义 ====
function Cb = calculateCb(A, d)
    Rf = 0.0738 * (d^1.021); 
    Zu = 9.95 * 0.8 * (Rf^0.882); 
    Cb = (A * Zu^5) / (1 + (A/0.3) * Zu^5);
end

function Sb = calculateSb(y, h1, a, Cb, P)
%     y(y == 0) = 0; 
    Sb = (((h1 - y) ./ y) .* (a / (h1 - a)) .^ (P)) * Cb;
end

