function [] = plot_perm_G
load('data/fit_result_NN.mat',...
    'rs','d_fitmodel_full','d_rs_add','d_rs_full','d_rs_multi');

rs(rs<0)=nan;d_rs_add(d_rs_add<0)=nan;d_rs_multi(d_rs_multi<0)=nan;d_rs_full(d_rs_full<0)=nan;

% startup
figure
set(gcf,'position',[300 200 1100 650]);
s(1) = subplot(2,3,1);
%CAESER_ARRAY M1
% R1
t = 1:size(rs,2);
% ADD
R = smooth(nanmean(d_rs_add(1:end,:),1))';
E = std(d_rs_add(1:end,:),1,'omitnan')/sqrt(size(d_rs_add(1:end,:),1));
plot(R,'linewidth',2,'color',[8 118 191]/255); hold on
fill([t fliplr(t)],[R+2*E fliplr(R-2*E)],[8 118 191]/255,'LineStyle','none');alpha(.3);
% COSINE
R = smooth(nanmean(rs(1:end,:),1))';
E = std(rs(1:end,:),1,'omitnan')/sqrt(size(rs(1:end,:),1));
plot(R,'linewidth',2,'color',[66 66 66]/255); hold on
fill([t fliplr(t)],[R+2*E fliplr(R-2*E)],[66 66 66]/255,'LineStyle','none');alpha(.3);
% FULL
R = smooth(nanmean(d_rs_full(1:end,:),1))';
E = std(d_rs_full(1:end,:),1,'omitnan')/sqrt(size(d_rs_full(1:end,:),1));
plot(R,'linewidth',2,'color',[217 83 25]/255); hold on
fill([t fliplr(t)],[R+2*E fliplr(R-2*E)],[217 83 25]/255,'LineStyle','none');alpha(.3);
% MULTI
R = smooth(nanmean(d_rs_multi(1:end,:),1))';
E = std(d_rs_multi(1:end,:),1,'omitnan')/sqrt(size(d_rs_multi(1:end,:),1));
plot(R,'linewidth',2,'color',[130 54 146]/255); hold on
fill([t fliplr(t)],[R+2*E fliplr(R-2*E)],[130 54 146]/255,'LineStyle','none');alpha(.3);

% effect size and significant
p=[];r=[];
for i=1:length(R)
multi = d_rs_multi(1:end,i);
add = d_rs_add(1:end,i); 
[p(i),~,stats] = signrank(add,multi,'method','approximate');
r(i) = stats.zval/sqrt(length(add));
end
% Msig = double(r<=-0.4); Msig(Msig<1)=nan;
% plot(Msig,'color',[130 54 146]/255,'linewidth',2);
% Asig = double(r>=0.4); Asig(Asig<1)=nan;
% plot(Asig,'color',[8 118 191]/255,'linewidth',2);
effect_plot(r)

box off
ax = gca;xlim([3.5,83.5]);ylim([0 1]);
set(gca,'xtick',3.5:10:83.5);
ax.XTickLabel = {'-0.8','','-0.4','','0','','0.4','','0.8'};
% ax.YTickLabel = {'','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};
title('Goodness of fit');ylabel('Adjusted R^2');xlabel('Time to MO(s)');

%Coeff
d_coff_full=[];
col=[21,112,177;29,153,29;255,119,0;118 113 113;141,194,211]/255;
for t=1:size(d_fitmodel_full,2)
    for n=1:size(d_fitmodel_full,1)
        d_coff_full(n,t,:) = d_fitmodel_full{n,t}.Coefficients.Estimate;
    end
end
a = squeeze(mean(abs(d_coff_full(1:end,:,[1 3 4 5])),1));
t = 1:size(rs,2);
coff = abs(d_coff_full(1:end,:,[1 3 4 5]));E=[];
for i=1:4
    E(:,i) = squeeze(std(coff(:,:,i),1)/sqrt(size(coff(:,:,i),1)));
end
s(3)=subplot(4,3,3);
plot(a(:,4),'color',col(4,:),'linewidth',2);box off;hold on
fill([t fliplr(t)],[a(:,4)'+2*E(:,4)' fliplr(a(:,4)'-2*E(:,4)')],col(4,:),'LineStyle','none');alpha(.3);
ax = gca;xlim([3.5,83.5]);title('Full model coefficient')
set(gca,'xtick',3.5:10:83.5);ylabel('Abs. Offset');ylim([0,0.5]);ax.YTick=[0,0.5];
% ax.XTickLabel = {'-0.8','','-0.4','','0','','0.4','','0.8'};
ax.XTick=[];
s(6)=subplot(4,3,6);
for i=1:3
plot(a(:,i),'linewidth',2,'color',col(i,:));hold on
fill([t fliplr(t)],[a(:,i)'+2*E(:,i)' fliplr(a(:,i)'-2*E(:,i)')],col(i,:),'LineStyle','none');alpha(.3);
end
box off
ax = gca;xlim([3.5,83.5]);ylim([0,0.3]);ax.YTick=[0,0.3];
set(gca,'xtick',3.5:10:83.5);ylabel('Abs. Coeff');xlabel('Time to MO(s)');
ax.XTickLabel = {'-0.8','','-0.4','','0','','0.4','','0.8'};
% s(6).Position(2)=0.58;

%Raster t=43
t_bin=43;
s(2)=subplot(2,3,2);

multi = d_rs_multi(1:end,t_bin);
add = d_rs_add(1:end,t_bin);

plot(add,multi,'.','markersize',10);
hold on
plot([0 1],[0,1],'r','linewidth',1.5);
axis([0 1 0 1]);
axis square
set(gca,'ytick',0:0.5:1);set(gca,'xtick',0:0.5:1);box off
xlabel('Additive');ylabel('Multiplicative');title('Model Preference');
annotation('textbox',[.42 .83 .2 .1],'String',['Wilcoxon signed rank test p=' char(vpa(p(t_bin),2))],'EdgeColor','none','FontSize',12);

%% plot perm after regress
load('data/Perm_result_fitnlmG1st.mat')
perm_rs(:,:,:,1) = d_rs_full;
perm_coff(:,:,:,:,1) = d_coff_full;
load('data/Perm_result_fitnlmGmul.mat')
perm_rs(:,:,:,2) = d_rs_full;
perm_coff(:,:,:,:,2) = d_coff_full;
load('data/Perm_result_fitnlmG2nd.mat')
perm_rs(:,:,:,3) = d_rs_full;
perm_coff(:,:,:,:,3) = d_coff_full;

clearvars -except perm_rs perm_coff t

rs_m = squeeze(nanmean(perm_rs,3));
coff_m = squeeze(mean(abs(perm_coff),4));


% M1
for p = 1:3
    rs = rs_m(1:end,:,p);
    coff = coff_m(1:end,:,:,p);
    
    %Coeff
    col=[21,112,177;29,153,29;255,119,0;118 113 113;141,194,211]/255;
    a = squeeze(mean(abs(coff(:,:,[1 3 4 5])),1));
    coff = abs(coff(:,:,[1 3 4 5]));E=[];
    for i=1:4
        E(:,i) = squeeze(std(coff(:,:,i),1)/sqrt(size(coff(:,:,i),1)));
    end
    
    s(6)=subplot(4,3,6);
    fill([t fliplr(t)],[a(:,p)'+2*E(:,p)' fliplr(a(:,p)'-2*E(:,p)')],col(p,:),'LineStyle','none');alpha(.5);
    hold on
    
end
s(6).Position(2)=0.58;
end