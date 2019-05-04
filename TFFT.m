function U = TFFT(x, y, d, Ui, k)
    [X, Y] = meshgrid(x, y);
    F0 = exp(1j*k*d)/(1j*d*(2*pi/k));
    fUi = fft2(Ui);
    F1 = exp(1j*k/2/d*(X.^2+Y.^2));
    fF1 = fft2(F1);
    Fuf = fUi.*fF1;
    U = F0.*fftshift(ifft2(Fuf));
end