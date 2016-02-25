function d = ddtf_dec_p(u, H)

N = size(H,1);
s = sqrt(N); % filter size

flag = mod(s,2)==0;

d = cell(N, 1);

u = u/s;
for i = 1:N
    
    h = H(:,i);
    h = reshape(h, s,s);
    if flag
        h = padarray(h, [1,1], 'post');
    end
    
    d{i} = imfilter(u, h, 'circular');
    
end
