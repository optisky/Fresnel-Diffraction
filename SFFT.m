function [x, y, U] = SFFT(xi, yi, d, Ui, k)
    % xi: x-coordinate of the object
    % yi: y-coordinate of the object
    % d:  distant of propagation
    % Ui: optical field the object
    % k:  wave number
    lam = 2*pi/k;
    Lx = length(xi)*lam*d/(xi(end)-xi(1));
    Ly = length(yi)*lam*d/(yi(end)-yi(1));
    x = linspace(-Lx, Lx, length(xi));
    y = linspace(-Ly, Ly, length(yi));
    [X, Y] = meshgrid(x, y);
    [Xi, Yi] = meshgrid(xi, yi);
    F0 = exp(1j*k*d)/(1j*lam*d)*exp(1j*k/2/d*(X.^2+Y.^2));
    F = Ui.*exp(1j*k/2/d*(Xi.^2+Yi.^2));
    Ff = fftshift(fft2(F));
    U = F0.*Ff;
end
