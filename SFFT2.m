function [x, y, U] = SFFT2(xi, yi, d, Ui, k, d2tod1)
    % Reference: 胡琪, 王喆, 刘洪顺, 薛智文, 邓家春. 基于单次傅里叶变换的分段衍射算法. 
    % 中国光学. 第11卷, 第4期, 2018.
    d1 = d / (1 + d2tod1);
    d2 = d - d1;
    [x1, y1, U1] = SFFT(xi, yi, d1, Ui, k);
    [x, y, U] = SFFT(x1, y1, d2, U1, k);
end