function u = ddtf_rec_p(d, H)

N = size(H,1);
s = sqrt(N); % filter size

flag = mod(s,2)==0;

u = 0;

for i = 1:N
    
    h = H(:,i);
    h = h(end:-1:1);
    h = reshape(h, s,s);
    if flag
        h = padarray(h, [1,1], 'pre');
    end
    
    u = u + imfilter(d{i}, h, 'circular');
    
end

u = u/s;