function [x] = pulsetrain(n, M, DC, offset)
    x = ones(1, length(n));
    for i=n
        if(mod(i+offset-1,M)>DC*M)
            x(i)=0;
        end
    end
end