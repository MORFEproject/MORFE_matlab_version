function Del = Delta(H)
Del = zeros(2*H+1);
for i = 1:H
    Del(2*i:2*i+1,2*i:2*i+1) = [0 -i; i 0];
end
end