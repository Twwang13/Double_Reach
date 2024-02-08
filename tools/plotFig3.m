function fig = plotFig3
load('data/exampleFit.mat');


%% Fig 3b
fig = temp(rs,rs_se,rs_sec);


%% plot Fig 3c
figure;
set(gcf,'position',[150 50 1300 842]);
n = 60; %example neuron 142
color_CY = [245, 124, 0; 162, 20, 47;183, 70, 255; 0, 114, 189; 119, 172, 48]/255;


add_model = d_fitmodel_add(n,:);z_add = [];
for t = 1:95
    fitmodel = add_model{t};
    z_add(t,:) = fitmodel(dr_am, dr_asm);%add
    f = fittype('a*cosd(x-b)+c');coff = coeffvalues(fitmodel);
    c = cfit(f,coff(1),coff(2),coff(3));
%     z_co_add(t,:) = feval(f,coff(1),coff(2),coff(3),co_am);
end 
mul_model = d_fitmodel_multi(n,:);z_mul=[];
for t = 1:95
    fitmodel = mul_model{t};
    z_mul(t,:) = fitmodel(dr_asm,dr_am);%multi
    f = fittype('a*cosd(x-b)+c');coff = coeffvalues(fitmodel);
    c = cfit(f,coff(1),coff(2),coff(3));
%     z_co_mul(t,:) = feval(f,coff(1),coff(2),coff(3),co_am);
end 
FULL_model = d_fitmodel_full(n,:);z_full=[];
for t = 1:95
    fitmodel = FULL_model{t};
    z_full(t,:) = fitmodel(dr_asm,dr_am);%FULL
    f = fittype('a*cosd(x-b)+c');coff = coeffvalues(fitmodel);
    c = cfit(f,coff(1),coff(2),coff(5));
%     z_co_full(t,:) = feval(f,coff(1),coff(2),coff(5),co_am);
end 
n_dsp = squeeze(dr_spm(n,:,:));z=n_dsp';
z_full(z_full<0)=0;z_add(z_add<0)=0;z_mul(z_mul<0)=0;% 将小于0的发放率设置为0
t = -0.85:0.02:1.03;%+0.15的修正
LW=2;
n_sp = squeeze(co_spm(n,:,:));z=[z,n_sp'];
for d=1:6
s = subplot(6,10,(d-1)*10+1);d=d-1;%first movement
xlim([-1,1]);ylim([-1,1]);axis square
hold on;
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'LineWidth',2,'LineStyle','--','EdgeColor',[200 200 200]/255);
scatter(0,0,100,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerFaceAlpha',.3);
scatter(cos(d*pi/3),sin(d*pi/3),100,'s','MarkerFacecolor',[0.3660 0.8740 0.5880],'MarkerEdgecolor',[0.4660 0.8740 0.5880]);
arrow([0,0],[cos(d*pi/3),sin(d*pi/3)],'width',2.3,'Length',10,'EdgeColor','k','FaceColor','k');
ax=gca;set(gca,'color','none');
ax.YAxis.Visible = 'off';   % 设置y轴不可见
ax.XAxis.Visible = 'off';   % 默认属性 on 表明可见
d=d+1;

s(d)=subplot(6,10,(d-1)*10+2);d=d-1;%second movements
xlim([-1,1]);ylim([-1,1]);axis square
hold on;
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'LineWidth',2,'LineStyle','--','EdgeColor',[200 200 200]/255);
scatter(0,0,100,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerFaceAlpha',.3);
arrow([cos(d*pi/3),sin(d*pi/3)],[cos((d+2)*pi/3),sin((d+2)*pi/3)],'width',1,'Length',10,'EdgeColor',color_CY(2,:),'FaceColor',color_CY(2,:));%
arrow([cos(d*pi/3),sin(d*pi/3)],[cos((d+4)*pi/3),sin((d+4)*pi/3)],'width',1,'Length',10,'EdgeColor',color_CY(4,:),'FaceColor',color_CY(4,:));%
arrow([cos(d*pi/3),sin(d*pi/3)],[cos((d+3)*pi/3),sin((d+3)*pi/3)],'width',1,'Length',10,'EdgeColor',color_CY(3,:),'FaceColor',color_CY(3,:));%
arrow([cos(d*pi/3),sin(d*pi/3)],[cos((d+1)*pi/3),sin((d+1)*pi/3)],'width',1,'Length',10,'EdgeColor',color_CY(1,:),'FaceColor',color_CY(1,:));%
arrow([cos(d*pi/3),sin(d*pi/3)],[cos((d-1)*pi/3),sin((d-1)*pi/3)],'width',1,'Length',10,'EdgeColor',color_CY(5,:),'FaceColor',color_CY(5,:));%
scatter(cos(d*pi/3),sin(d*pi/3),100,'s','MarkerFacecolor',[0.3660 0.8740 0.5880],'MarkerEdgecolor',[0.4660 0.8740 0.5880]);
ax=gca;set(gca,'color','none');
ax.YAxis.Visible = 'off';   % 设置y轴不可见
ax.XAxis.Visible = 'off';   % 默认属性 on 表明可见
d=d+1;

%RAW
s(d+6) = subplot(6,5,(d-1)*5+2);
plot(t,z(:,d)','color',color_CY(2,:),'LineWidth',LW);hold on
plot(t,z(:,d+6)','color',color_CY(4,:),'LineWidth',LW);
plot(t,z(:,d+12)','color',color_CY(3,:),'LineWidth',LW);
plot(t,z(:,d+18)','color',color_CY(1,:),'LineWidth',LW);
plot(t,z(:,d+24)','color',color_CY(5,:),'LineWidth',LW);box off
% plot(t,z(:,d+30)','k','LineWidth',LW);
ylim([0 1]);xlim([-0.8,0.6]);%ax=gca;ax.XAxis.Visible = 'off'; 
%ADD
s(d+12) = subplot(6,5,(d-1)*5+3);
plot(t,z_add(:,d)','color',color_CY(2,:),'LineWidth',LW);hold on
plot(t,z_add(:,d+6)','color',color_CY(4,:),'LineWidth',LW);
plot(t,z_add(:,d+12)','color',color_CY(3,:),'LineWidth',LW);
plot(t,z_add(:,d+18)','color',color_CY(1,:),'LineWidth',LW);
plot(t,z_add(:,d+24)','color',color_CY(5,:),'LineWidth',LW);box off
% plot(t,z_co_add(:,d)','k','LineWidth',LW);
ylim([0 1]);xlim([-0.8,0.6]);ax=gca;ax.YAxis.Visible = 'off';   % 设置y轴不可见
% ax.XAxis.Visible = 'off'; 
%MUL
s(d+18) = subplot(6,5,(d-1)*5+4);
plot(t,z_mul(:,d)','color',color_CY(2,:),'LineWidth',LW);hold on
plot(t,z_mul(:,d+6)','color',color_CY(4,:),'LineWidth',LW);
plot(t,z_mul(:,d+12)','color',color_CY(3,:),'LineWidth',LW);
plot(t,z_mul(:,d+18)','color',color_CY(1,:),'LineWidth',LW);
plot(t,z_mul(:,d+24)','color',color_CY(5,:),'LineWidth',LW);box off
% plot(t,z_co_mul(:,d)','k','LineWidth',LW);
ylim([0 1]);xlim([-0.8,0.6]);ax=gca;ax.YAxis.Visible = 'off';   % 设置y轴不可见
% ax.XAxis.Visible = 'off'; 
%FULL
s(d+24) = subplot(6,5,(d-1)*5+5);
plot(t,z_full(:,d)','color',color_CY(2,:),'LineWidth',LW);hold on
plot(t,z_full(:,d+6)','color',color_CY(4,:),'LineWidth',LW);
plot(t,z_full(:,d+12)','color',color_CY(3,:),'LineWidth',LW);
plot(t,z_full(:,d+18)','color',color_CY(1,:),'LineWidth',LW);
plot(t,z_full(:,d+24)','color',color_CY(5,:),'LineWidth',LW);box off
% plot(t,z_co_full(:,d)','k','LineWidth',LW);
ylim([0 1]);xlim([-0.8,0.6]);ax=gca;ax.YAxis.Visible = 'off';   % 设置y轴不可见
% ax.XAxis.Visible = 'off'; 

end
subplot(6,5,2);title('PSTHs');subplot(6,5,3);title('Additive');subplot(6,5,4);title('Multiplicative');subplot(6,5,5);title('FULL');
subplot(6,5,27);xlabel('Time to MO (s)');subplot(6,5,28);xlabel('Time to MO (s)');subplot(6,5,29);xlabel('Time to MO (s)');subplot(6,5,30);xlabel('Time to MO (s)');
% subplot(6,5,12);ylabel('Normalized FR','FontWeight','bold','position',[-1.3 -0.2]);
subplot(6,10,1);title('1st Reach');subplot(6,10,2);title('2nd Reach');
end