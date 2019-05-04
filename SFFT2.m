function [x, y, U] = SFFT2(xi, yi, d, Ui, k, d2tod1)
    % Reference: ����, ����, ����˳, Ѧ����, �˼Ҵ�. ���ڵ��θ���Ҷ�任�ķֶ������㷨. 
    % �й���ѧ. ��11��, ��4��, 2018.
    d1 = d / (1 + d2tod1);
    d2 = d - d1;
    [x1, y1, U1] = SFFT(xi, yi, d1, Ui, k);
    [x, y, U] = SFFT(x1, y1, d2, U1, k);
end