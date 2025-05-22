clc;
clear all;
close all;

% 生成单宽流量q
% 生成时间序列
time = linspace(0, 24, 1441); % 0到24小时，1分钟为间隔

% 生成流量数据
peak_time = 4; % 峰值时间为第4小时
peak_flow = 1000; % 峰值流量

% 高斯分布峰值
flow_peak = peak_flow * exp(-0.5 * ((time - peak_time) / 1.4).^2);

% 平滑连接峰值后的下降
decay_indices = time > peak_time & time <= 24*60;
decay_time = time(decay_indices);
flow_decay = peak_flow * exp(-0.1 * ((decay_time - peak_time) / 2).^2);

% 合并流量数据
flow = zeros(size(time));
flow(1:numel(flow_peak)) = flow_peak;
flow(decay_indices) = flow_decay;

% 计算剪切流速
% 参数
s = 0.0015; % 根据雅江中游（加查-朗县）的河道坡度，设定河道坡度为0.0015
g = 9.81; % 重力加速度
u_shear = (flow * g * s / 8.3 ).^(1/3);
v = 0.000001; 

% 计算bed shear stress
rho = 1000; % water density
rhos = 2650;
R = (rhos - rho) / rho;
tau = rho * (u_shear.^2); % tau is shear stress

% 计算平均速度
% 设定(8/fc)^0.5=8.3
U = 8.3 * u_shear;

% 创建图形窗口
figure;
set(gcf, 'Position', [100, 100, 1200, 900]);

% 子图1：绘制流量曲线
subplot(3, 2, 1); % 2行2列的子图，当前绘制第1个子图
plot(time, flow,'Color', [62/255 38/255 169/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Discharge per unit width {\itq} (m^2/s)', 'FontSize', 14);
xlim([0 24]);
ylim([0 1200]);
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])
%grid on; % 添加网格

% 添加辅助线
line([4 4], [0 1000], 'linestyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
line([0 4], [1000 1000], 'linestyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

% 子图2：绘制平均流速
subplot(3, 2, 2); % 2行2列的子图，当前绘制第1个子图
plot(time, U,'Color', [244/255 120/255 58/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Mean velocity {\itu} (m/s)', 'FontSize', 14);
xlim([0 24]);
ylim([0 15]); % 自动调整 y 轴范围
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])
% grid on; % 添加网格
% 在左上角添加序号 "(b)"
%text(0.5, max(U) * 1.2, '(b)', 'FontSize', 16, 'FontWeight', 'bold');

% 子图3：绘制剪切流速
subplot(3, 2, 3); % 2行2列的子图，当前绘制第3个子图
plot(time, u_shear,'Color', [127/255 127/255 127/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Shear velocity {\itu_*} (m/s)', 'FontSize', 14);
xlim([0 24]);
ylim([0 1.5]); 
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])
ytickformat('%.1f');
% grid on; % 添加网格
% 在左上角添加序号 "(c)"
%text(0.5, max(u_shear) * 1.2, '(c)', 'FontSize', 16, 'FontWeight', 'bold');


% 子图4：绘制剪切力曲线
subplot(3, 2, 4); % 2行2列的子图，当前绘制第4个子图
plot(time, tau,'Color', [12/255 172/255 232/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Bed shear stress {\it\tau} (N/m^2)', 'FontSize', 14);
xlim([0 24]);
ylim([0 2000]); % 自动调整 y 轴范围
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])
% grid on; % 添加网格
% 在左上角添加序号 "(d)"
%text(0.5, max(tau) * 1.3, '(d)', 'FontSize', 16, 'FontWeight', 'bold');

% % 子图5：绘制水流雷诺数
% % 计算水深
% h = (u_shear .^2) / g / s;
% % 计算水流Renolds Number
% Re = U .* h ./ v;
% subplot(3, 2, 5); % 2行2列的子图，当前绘制第4个子图
% plot(time, Re,'Color', [69/255 141/255 252/255], 'LineWidth', 3.5); 
% xlabel('Time {\itT} (hr)', 'FontSize', 14);
% ylabel('Flow Renolds number {\itRe} (-)', 'FontSize', 14);
% xlim([0 24]);
% ylim([1e2 1e10]); 
% set(gca, 'YScale', 'log'); 
% set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
% set(gca, 'TickLength', [0.025 0.025]);
% 
% 
% % 子图6：绘制水流弗劳德数
% % 计算Froude Number
% % W = 2000; % channel width
% % h = (u_shear .^2) / g / s;
% % R = h .* W ./ (2 * h + W); % hydraulic radius
% % Fr = U ./ sqrt(g .* R); 
% 
% h = (u_shear .^2) / g / s;
% Fr = U ./ sqrt(g .* h); 
% subplot(3, 2, 6); 
% plot(time, Fr,'Color', [196/255 147/255 255/255], 'LineWidth', 3.5); 
% xlabel('Time {\itT} (hr)', 'FontSize', 14);
% ylabel('Froude Number {\itFr} (-)', 'FontSize', 14);
% xlim([0 24]);
% ylim([0 1]); % 自动调整 y 轴范围
% set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
% set(gca, 'TickLength', [0.02 0.025])
% 
% export_fig HydraulicVariable.jpg -r1000 -transparent;




% 计算水深
h = (u_shear .^2) / g / s;
figure;

% 水深图
plot(time, h, 'LineWidth', 4, 'Color', [69/255 141/255 252/255]); % 灰色曲线
xlabel('Time {\itT} (hr)', 'FontSize', 28);
ylabel('Water depth {\ith} (m)', 'FontSize', 28);
xlim([0 24]);
ylim([0 120]); % 自动调整 y 轴范围
set(gca, 'xtick', 0:2:24, 'FontSize', 20); 
set(gca, 'TickLength', [0.01 0.01]);
grid on; % 添加网格
hold on;

% 计算特定时间点
specific_times = [120, 240, 360, 480, 600, 720, 840, 960]; % 分钟
specific_times_hours = specific_times / 60; % 转换为小时

% 计算对应时间点的水深值和平均流速值
h_specific_waterdepth = interp1(time, h, specific_times_hours, 'linear');
U_specific_averagevelocity = interp1(time, U, specific_times_hours, 'linear');


% 输出水深和平均流速
disp('特定时间点的水深值 (m):');
disp('时间 (小时)    水深 h (m)    平均流速 U (m/s)');        
disp([specific_times_hours', h_specific_waterdepth', U_specific_averagevelocity']);
          
% 自定义点的颜色
colors = [
    063/255, 007/255, 077/255; 
    068/255, 057/255, 131/255;
    047/255, 105/255, 140/255; 
    032/255, 145/255, 141/255;
    054/255, 184/255, 120/255; 
    108/255, 205/255, 90/255; 
    143/255, 215/255, 067/255; 
    242/255, 226/255, 048/255; 
];

% 在水深图上绘制垂直辅助线和不同颜色的交点
for i = 1:length(specific_times_hours)
    % 绘制垂直辅助线
    line([specific_times_hours(i), specific_times_hours(i)], [0, h_specific_waterdepth(i)], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    % 绘制不同颜色的交点
    plot(specific_times_hours(i), h_specific_waterdepth(i), 'o', 'MarkerSize', 18, 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k'); 
end
hold off;

% 调整图形窗口大小
set(gcf, 'Position', [100, 100, 1500, 400]);
%export_fig WaterDepth.tif -r1000 -transparent;


% % 绘制颗粒雷诺数
% figure;
% set(gcf, 'Position', [100, 100, 800, 350]);
% 
% % 生成 viridis 调色板
% map1 = flipud(mymap('viridis'));  % 确保不翻转 colormap
% colormap(gca, map1);  % 让 colormap 仅影响当前子图
% 
% % 生成颜色索引
% idx = linspace(1, size(map1, 1), length(D1) + length(D2));
% colors = map1(round(idx), :);
% 
% % 绘制 D1 对应的颗粒雷诺数曲线
% hold on;
% for i = 1:length(D1)
%     plot(time, Rep1(:, i), 'LineWidth', 3, 'Color', colors(i, :));
% end
% 
% % 绘制 D2 对应的颗粒雷诺数曲线
% for i = 1:length(D2)
%     plot(time, Rep2(:, i), 'LineWidth', 3, 'Color', colors(length(D1) + i, :));
% end
% 
% % 设置子图属性
% xlabel('Time (h)', 'FontSize', 14);
% ylabel('Particle Reynolds number (-)', 'FontSize', 14);
% xlim([0 24]);
% ylim([0 max([Rep1(:); Rep2(:)]) * 1.1]); 
% set(gca, 'xtick', 0:4:24, 'FontSize', 14, 'TickLength', [0.02 0.025]);
% grid on; box on;
%text(0.04, 0.95, '(a)', 'Units', 'Normalized', 'FontSize', 16, 'FontWeight', 'bold');
% 
% % 添加 colorbar 并正确设置刻度
% cbar = colorbar;
% cbar.Label.String = 'Particle Size (mm)'; 
% cbar.FontSize = 14;
% D_all = [D1, D2];
% caxis([min(D_all), max(D_all)]); % 设定 colorbar 范围
% cbar.Ticks = linspace(min(D_all), max(D_all), numel(D_all)); % 自定义刻度位置
% cbar.TickLabels = arrayfun(@(x) num2str(x * 1000), D_all, 'UniformOutput', false); % 自定义刻度标签
% cbar.TickLength = 0.03;
% hold off;





