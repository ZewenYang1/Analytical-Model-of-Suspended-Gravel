clc;
clear all;
close all;

% 基本常量
g = 9.81; % 重力加速度
K = 0.41; % 卡门常量
rho = 1000; % 水密度
rhos = 2650; % 颗粒密度
R = (rhos - rho) / rho; 
v = 1e-6; % 水的运动粘度
Mu = 0.001;

% 颗粒粒径范围 (m)
D = [0.001, 0.005, 0.010, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05];

% 8 组水深 (h) 和平均流速 (U)

h_values = [50.42, 99.55, 93.12, 76.25, 54.63, 34.26, 18.80, 9.03];
U_values = [7.15, 10.05, 9.72, 8.79, 7.44, 5.89, 4.37, 3.02];

% 配色方案
colors = parula(length(D)); 
colors = flipud(colors); % 反转颜色顺序


figure;
map = flipud(mymap('viridis'));
colormap(gca, map);
colormap(flipud(parula)); % 确保 colormap 与曲线颜色一致
set(gcf, 'Position', [100, 100, 1400, 1000]);

for i = 1:8
    h = h_values(i);
    U = U_values(i);
    
    % 垂直剖面位置
    y = linspace(0.05, h, 100);
    y_h = y / h;
    
    % 计算剪切流速
    u_shear = U / 8.3;
    
    % 计算湍动能耗散率
    Epsilon = ((9.8 * exp(-3 * y / h)) ./ ((y / h) .^ 0.5)) .* ((u_shear ^ 3) / h);
    
    % 计算 y1
    y1 = y * u_shear / v;
    
    % 计算方均根流速
    u_rootmean = 2.3 * u_shear .* exp(-y / h) .* (1 - exp(-y1 / 10)) + 0.3 * y1 .* exp(-y1 / 10);
    
    % 计算 Taylor length scale
    Lambda = sqrt(15 * v * (u_rootmean.^2) ./ Epsilon);
    
    % 创建子图
    subplot(4,2,i);
    set(gca, 'FontSize', 16);
    hold on;
    
    % 遍历不同粒径 D
    for j = 1:length(D)
        % 计算沉积流速
        Ws = 0.51 * (v / D(j)) * ((D(j)^3) * g * R / (v^2))^0.553; % Carling et al., 2020
        
        % 计算颗粒响应时间
        Tp = Ws / g;
        
        % 计算湍流时间尺度
        Tk = Lambda ./ u_rootmean;
        
        % 计算斯托克斯数
        St = Tp ./ Tk;
        
        % 绘制曲线
        plot(St, y_h, 'Color', colors(j, :), 'LineWidth', 2);
    end
    
    % 在子图上添加 St=1 和 St=9 的垂直线
    plot([1, 1], [0, 1], 'g--', 'LineWidth', 2); % St = 1 
    plot([9, 9], [0, 1],'--','Color', [0.5, 0.5, 0.5], 'LineWidth', 2); % St = 9 
    
    % 轴标签
    xlabel('Stokes number {\itS_t} (-)', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k'); 
    ylabel('{\ity/h} (-)', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
      
    % 设定轴的方向
    set(gca, 'YDir', 'normal', 'TickLength', [0.02 0.02]);
    set(gca, 'FontSize', 20);
    grid on;
    xlim([0 40]);
    ylim([0 1]);
    hold off;
    box on;
end

% 添加 colorbar 并映射 D 值
c = colorbar('eastoutside');
caxis([min(D) max(D)]);  
c.Label.String = 'Grain size {\itD} (mm)';  
c.Label.FontSize = 20;
c.Label.FontName = 'Arial';

% 设置 colorbar 刻度
tick_positions = linspace(min(D), max(D), length(D));
c.Ticks = tick_positions; 
c.TickLabels = arrayfun(@(x) num2str(round(x * 1000)), D, 'UniformOutput', false); % 显示为 mm
c.Position = [0.92, 0.121, 0.018, 0.79];
c.FontSize = 20;
c.Position = [c.Position(1), c.Position(2) + 0.035, c.Position(3), c.Position(4)-0.02];
%export_fig StProfileDistribution.jpg -r1000 -transparent;


% 额外的 figure 用于绘制 u' 在水深剖面上的分布
figure;
set(gcf, 'Position', [100, 100, 1400, 1000]);
for i = 1:8
    h = h_values(i);
    
    % 计算 u_rms 的分布
    y = linspace(0.05, h, 100);
    y_h = y / h;
    u_shear = U_values(i) / 8.3;
    y1 = y * u_shear / v;
    u_rootmean = 2.3 * u_shear .* exp(-y / h) .* (1 - exp(-y1 / 10)) + 0.3 * y1 .* exp(-y1 / 10);
    
    % 创建子图
    subplot(4,2,i);
    plot(u_rootmean, y_h, 'b', 'LineWidth', 4,'Color', [112/255, 48/255, 160/255]);
    set(gca, 'FontSize', 20);
    xlabel("u_{rms} (m/s)", 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    ylabel('y/h', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    %title(['h = ', num2str(h), 'm'], 'FontSize', 16, 'FontName', 'Arial', 'Color', 'k');
    set(gca, 'YDir', 'normal');
    grid on;
end
