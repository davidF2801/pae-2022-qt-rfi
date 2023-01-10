function [x_filt] = lowpass_FFT(x, fcut, fs, Nfft)
    filter = zeros(1, Nfft);
    cut = fcut*Nfft/fs;
    filter(1:cut) = ones(1,cut);
    filter((Nfft-cut+1):Nfft) = ones(1,cut);
    fft_x = fft(x,Nfft);
    fft_x_filt = fft_x.*filter;
    x_filt = ifft(fft_x_filt, Nfft);
end