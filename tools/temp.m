function fig = temp(rs,rs_se,rs_sec)
fig = figure
t = -0.85:0.02:1.03;
R = mean(rs,1);
E = std(rs,1,'omitnan')/sqrt(size(rs,1));
plot(t,R,'linewidth',2,'color',[66/255 66/255 66/255]);hold on
fill([t fliplr(t)],[R+2*E fliplr(R-2*E)],[66/255 66/255 66/255],'LineStyle','none');alpha(.3);


R = mean(rs_se,1);
E = std(rs_se,1,'omitnan')/sqrt(size(rs_se,1));
plot(t,R,'linewidth',2,'color',[201/255 92/255 46/255]);hold on
fill([t fliplr(t)],[R+2*E fliplr(R-2*E)],[201/255 92/255 46/255],'LineStyle','none');alpha(.3);

R = mean(rs_sec,1);
E = std(rs_sec,1,'omitnan')/sqrt(size(rs_sec,1));
plot(t,R,'linewidth',2,'color',[48/255 112/255 183/255]);hold on
fill([t fliplr(t)],[R+2*E fliplr(R-2*E)],[48/255 112/255 183/255],'LineStyle','none');alpha(.3);

ax = gca;xlim([-0.8,0.8]);ylim([0 1]);axis square
set(gca,'xtick',-0.8:0.2:0.8);
ax.XTickLabel = {'-0.8','-0.6','-0.4','-0.2','MO','0.2','0.4','0.6','0.8'};
box off
% ax.YTickLabel = {'','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};
xlabel('time (s)');ylabel('R^2');title('Goodness of fit');

end