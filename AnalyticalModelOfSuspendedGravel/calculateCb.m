% 计算参考浓度的函数
function Cb = calculateCb(d)
    A = 1.3e-7;
    Rf = 0.0738 * (d.^1.021);
    Zu = 9.95 * 0.8 * (Rf.^0.882);
    Cb = (A * Zu.^5) ./ (1 + (A/0.3) * Zu.^5);
end