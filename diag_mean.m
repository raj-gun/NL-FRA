function out=diag_mean(mat)
len_diag=length(diag(mat));
out=zeros(len_diag,1);
for i=0:len_diag-1
    out(i+1) = mean(diag(mat,i));
end
end