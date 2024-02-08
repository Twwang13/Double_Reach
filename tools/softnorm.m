function [xx] = softnorm(data)
tmp = reshape(data,[size(data,1),size(data,2)*size(data,3)]);%NKT
tmp = tmp./repmat(range(tmp')'+5,[1 size(tmp,2)]);
xx = reshape(tmp,[size(data,1),size(data,2),size(data,3)]);

end