function index = find_index_VinA(V, A)
AminusV_norm = sum(abs(A-V),1);
index = find(AminusV_norm == 0);
end