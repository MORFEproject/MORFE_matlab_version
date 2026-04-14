function A = IVec(Vec_A, n)
% 将矢量化的Vec_A转化为原矩阵A, 输入量还需要知道A的每个列向量的元素数目
A = reshape(Vec_A,n,[]);
end