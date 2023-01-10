function qsignal = quantize(signal, N, Vmin, Vmax)
    %N is the number of quantization bits of the ADC
    %A is the input value corresponding to the maximum output level
    %signal is the "not quantized" input signal
    %qsignal is the quantized output signal
    num_levels = (Vmax-Vmin)/(2^N);
    qsignal = round(signal/num_levels)*num_levels;
end