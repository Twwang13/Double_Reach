function  [direction_ave] = direction_ave_cal(xPoints, yPoints, binsize, firstBin, lastBin, mBin, nBin)
edges =  firstBin : binsize : lastBin;	
ntrials = length(unique(yPoints));
h = (histcounts(xPoints,edges));
x = h*1/binsize/ntrials;
center = abs(firstBin)/binsize;
s_x = smoothdata(x, 'gaussian', 5);
direction_ave = mean(s_x(center+mBin/binsize:center+nBin/binsize));
end