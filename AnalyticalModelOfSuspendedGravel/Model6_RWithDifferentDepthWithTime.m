% 该代码用于判断悬移质水流中是否会出现前锋（Front）（参考自Bennett et al.,2014 图12）
clc;
clear all;
close all force;

% 1.时间与流量构建 
time = linspace(0, 24, 1441);  % 每分钟一个点
peak_time = 4;     % 峰值流量发生时间（小时）
peak_flow = 1000;  % 峰值流量（m2/s）

flow_peak = peak_flow * exp(-0.5 * ((time - peak_time) / 1.4).^2);  % 高斯分布上升段
decay_indices = time > peak_time & time <= 24;
decay_time = time(decay_indices);
flow_decay = peak_flow * exp(-0.1 * ((decay_time - peak_time) / 2.0).^2);  % 衰减段

flow = zeros(size(time));
flow(1:numel(flow_peak)) = flow_peak;
flow(decay_indices) = flow_decay;

% 2.常量定义
s = 0.0015;
rho = 1000;     % 水密度
rhos = 2650;    % 颗粒密度
R = (rhos - rho) / rho;
g = 9.81;       % 重力加速度
v = 1e-6;       % 动力粘度
karmanconstant = 0.41;

% 3.粒径定义
D1 = 0.001:0.001:0.005;
D2 = 0.01:0.01:0.05;
D_values = [D1, D2];  % 合并所有粒径
nD = length(D_values);

% 4.剪切流速
u_shear = (flow * g * s / 8.3 ).^(1/3);  % 剪切流速
tau = rho * u_shear.^2;

% 5.建立颜色映射
cmap = parula(nD); 
colormap(cmap);

% 6.设定不同无量纲高度a/h 值
k_list = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];

% 7.绘图
figure;
set(gcf, 'Position', [200, 200, 1200, 950]); 

for idx = 1:length(k_list)
    k_ratio = k_list(idx);

    subplot(3, 2, idx); hold on; box on;
    %title(['a/h = ', num2str(k_ratio)], 'FontSize', 16);

    for i = 1:nD
        D = D_values(i); 

        % 沉降流速
        Ws = 0.51 * (v / D) * ((D^3 * g * R) / v^2)^0.553;

        % 水深
        h = (u_shear.^2) / g / s;

        % 参考高度
        a = k_ratio * h;

        % z+ 参数
        z = a .* u_shear / v;

        % 均方根速度
        u_rootmean = u_shear .* 2.3 .* exp(-a./h) .* (1 - exp(-z ./ 10)) + 0.3 .* z .* exp(-z ./ 10);

        % ==== R_star 计算 ====
        R_star = Ws ./ u_rootmean;

        % 绘图
        plot(time, R_star, 'Color', cmap(i,:), 'LineWidth', 2.5);
    end
    xlabel('Time {\itT} (hr)', 'FontSize', 18);
    ylabel('{\itR} (-)', 'FontSize', 18);
    xlim([0, 24]);
    ylim([0, 60]);
    xticks(0:4:24);
    set(gca, 'FontSize', 18);
    set(gca, 'TickLength', [0.015, 0.015]);
    grid on;
end

hcb = colorbar('Position', [0.93 0.133 0.0185 0.791]);
caxis([1, nD]);
hcb.Ticks = linspace(1, nD, nD);
hcb.TickLabels = arrayfun(@(x) sprintf('%.d', x*1000), D_values, 'UniformOutput', false);
ylabel(hcb, 'Grain size {\itD} (mm)', 'FontSize', 18);
hold off;

%export_fig R(-).jpg -r1000 -transparent;
