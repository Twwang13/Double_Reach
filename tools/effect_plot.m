function effect_plot(r)

im = imagesc(r);

Cmin = abs(fix(min(r)*100))-20;Cmax = abs(fix(max(r)*100))-20;
tmp = [130 54 146;8 118 191]/255;
mmp = [linspace(tmp(2,1),1,Cmax)',linspace(tmp(2,2),1,Cmax)',linspace(tmp(2,3),1,Cmax)'];
mmp = [mmp;[linspace(1,1,40)',linspace(1,1,40)',linspace(1,1,40)']];
mmp = [mmp;[linspace(1,tmp(1,1),Cmin)',linspace(1,tmp(1,2),Cmin)',linspace(1,tmp(1,3),Cmin)']];
colormap(flipud(mmp))
im.YData = 1.45;
end