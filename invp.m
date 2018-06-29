function ip = invp(p)

ip = zeros(length(p),1);
for i = 1:length(p)
    ip(p(i)) = i;
end