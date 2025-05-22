% 该代码用于判断悬移质水流中是否会出现前锋（Front）（参考自Bennett et al.,2014 图12）

clc;
clear all;
close all;

% ==== 1. 时间与流量构建 ====
time = linspace(0, 24, 1441);  % 每分钟一个点
peak_time = 4;     % 峰值流量发生时间（小时）
peak_flow = 1000;  % 峰值流量（m?/s）

flow_peak = peak_flow * exp(-0.5 * ((time - peak_time) / 1.4).^2);  % 高斯分布上升段
decay_indices = time > peak_time & time <= 24;
decay_time = time(decay_indices);
flow_decay = peak_flow * exp(-0.1 * ((decay_time - peak_time) / 2).^2);  % 衰减段

flow = zeros(size(time));
flow(1:numel(flow_peak)) = flow_peak;
flow(decay_indices) = flow_decay;

% ==== 2. 常量定义 ====
s = 0.0015;
rho = 1000;     % 水密度
rhos = 2650;    % 颗粒密度
R = (rhos - rho) / rho;
g = 9.81;       % 重力加速度
v = 1e-6;       % 动力粘度
karmanconstant = 0.41;

% ==== 3. 粒径定义 ====
D1 = 0.001:0.001:0.005;
D2 = 0.01:0.01:0.05;
D_all = [D1, D2];  % 合并所有粒径
nD = length(D_all);

% ==== 4. 剪切流速 ====
u_shear = (flow * g * s / 8.3 ).^(1/3);  % 剪切流速
tau = rho * u_shear.^2;

% ==== 5. 建立颜色映射 ====
%colors = parula(length(D)); 
map = flipud(parula(length(D_all))); 
color_all = map(round(linspace(1, size(map, 1), nD)), :);
% 生成颜色索引
idx1 = linspace(1, size(map, 1) / 2, length(D1)); % D1 使用 colormap 前半部分
idx2 = linspace(size(map, 1) / 2 + 1, size(map, 1), length(D2)); % D2 使用 colormap 后半部分
colors1 = map(round(idx1), :);
colors2 = map(round(idx2), :);


% ==== 6. 绘图 ====
% ==== a/h 比例列表 ====
k_list = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];

% ==== 绘图 ====
figure;
set(gcf, 'Position', [200, 150, 1200, 950]); 

% ==== 遍历每个 a/h 比例 ====
for idx = 1:length(k_list)
    k_ratio = k_list(idx);  % 当前的 a/h 比例

    % ==== 创建对应子图 ====
    subplot(3, 2, idx); 
    hold on; box on;
    %title(['a/h = ', num2str(k_ratio)], 'FontSize', 16);

    % ==== 遍历每个粒径值 ====
    for i = 1:nD
        D = D_all(i);  % 当前粒径 (m)

        % ==== 计算沉降速度 Ws（m/s），公式来自 Soulsby (1997) ====
        Ws = 0.51 * (v / D) * ((D^3 * g * R) / v^2)^0.553;

        % ==== 水深 h（m），由剪切流速反推得出 ====
        h = (u_shear.^2) / g / s;

        % ==== 参考高度 a（m） ====
        a = k_ratio * h;

        % ==== z+ 参数（无量纲） ====
        z = a .* u_shear / v;

        % ==== 计算均方根流速 u_rootmean（m/s） ====
        u_rootmean = u_shear .* 2.3 .* exp(-a./h) .* ...
                     (1 - exp(-z ./ 10)) + 0.3 .* z .* exp(-z ./ 10);

        % ==== 湍动能耗散率 ε（m?/s?） ====
        Epsilon = (9.8 * exp(-3 * a ./ h) .* (u_shear.^3)) ./ ...
                  (((a ./ h) .^0.5) .* h);

        % ==== Taylor 尺度 λ（m） ====
        Lambda = sqrt(15 * v * (u_rootmean.^2) ./ Epsilon);

        % ==== 无量纲粒径 d，用于参考浓度计算 ====
        d = D .* ((R * g).^(1/3)) ./ (v^(2/3));

        % ==== 参考浓度 Cb（-） ====
        Cb = calculateCb(d);  % 调用外部函数

        % ==== R_star，无量纲沉降速率 ====
        R_star = Ws ./ u_rootmean;

        % ==== G_star，无量纲浓度动能参数 ====
        G_star = g * R * Lambda .* Cb ./ (u_rootmean.^2);

        % ==== G_star / R_star^2（目标比值） ====
        ratio1 = G_star ./ (R_star .^2);

        % ==== 绘图，每条曲线表示一个粒径 ====
        plot(time, ratio1, 'Color', color_all(i,:), 'LineWidth',2.5);
    end

    % ==== 坐标轴与格式设置 ====
    xlabel('Time {\itT} (hr)', 'FontSize', 18);
    ylabel('{\itG/R^2} (-)', 'FontSize', 18);
    xlim([0, 24]);
    ylim([0, 6]);  % 可按结果调整
    xticks(0:4:24);
    set(gca, 'FontSize', 18);
    set(gca, 'TickLength', [0.015, 0.015]);
    grid on;

    % ==== 添加三条参考线（对应前锋形成判别） ====
    plot([0, 24], [2.2, 2.2], 'g--', 'LineWidth', 1.5);
    plot([0, 24], [0.5, 0.5], 'Color', [0.5 0.5 0.5],'LineStyle', '--', 'LineWidth', 1.5);                  
    plot([0, 24], [0.032, 0.032], 'y--', 'LineWidth', 1.5);    
end


% ==== 添加 colorbar ====
colormap(flipud(parula)); 
c = colorbar('eastoutside');
caxis([min(D_all) max(D_all)]);

c.Label.String = 'Grain size {\itD} (mm)';  
c.Label.FontSize = 18;
c.Label.FontName = 'Arial';

% colorbar 刻度设置
tick_positions = linspace(min(D_all), max(D_all), length(D_all));
c.Ticks = tick_positions;
c.TickLabels = arrayfun(@(x) num2str(round(x * 1000)), D_all, 'UniformOutput', false); 
c.Position = [0.93 0.122 0.0187 0.803];

%export_fig GandR^2.jpg -r1000 -transparent;

function Cb = calculateCb(d)
    A = 1.3e-7;
    Rf = 0.0738 * (d.^1.021);
    Zu = 9.95 * 0.8 * (Rf.^0.882);
    Cb = (A * Zu.^5) ./ (1 + (A/0.3) * Zu.^5);
end
