% Generates a burst of pulses with period M, duty cycle DC and offset, where n is an array


function [x] = pulsetrain(n, M, DC, offset)
    %M      : lenght of a period
    %DC     : duty cycle (0<=DC<=1)
    %offset : number of points before the pulse
    x = ones(1, length(n));
    for i=n
        if(mod(i+offset-1,M)>=DC*M)
            x(i)=0;
        end
    end
end