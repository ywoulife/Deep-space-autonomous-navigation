%产生均值为0 的高斯白噪声
function out=wgn(row,column,~,~)
out=randn(row,column);
