clear all
L = 12; % time constant mumean resolution
S = 40000; % initial number of x samples
s = (1:S); % sample index update)
A0 = 1;
fs = 5e6;
beta = sqrt(30);
k = zeros(1,1000);
for i=1:10000
    noise = A0*normrnd(0,1,[1,S]); 
    r1 = noise;
    p1 = 10*log10(abs(r1).^2);
    k(i) = std(p1)/mad(p1);
end
c = mean(k);




