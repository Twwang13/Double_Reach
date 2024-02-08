function M = PCALDA_NC

load('data/NKT_158.mat');
clearvars -except NKT event name

bin_w = 300; 
stp = 300; 

GM = fix(1000*(event(:,4)-event(:,3)));
Atime = [-1300:1000]+GM;
X = [];
for i = 1:size(NKT,2)
    p_point = [find(Atime(i,:)==0)-599,find(Atime(i,:)==0)];
    X(:,i,:) = NKT(:,i,p_point(1):p_point(2));
end
bin = [];
for i =1:fix((size(X,3)-bin_w)/stp+1)
    bin(:,:,i)=sum(X(:,:,stp*(i-1)+1:stp*(i-1)+bin_w),3);
end
t = (-600:stp:0-bin_w)+(bin_w/2);

data=bin;
DATA = data(41:end,:,:);

J_data = [];
condition = [];
conv_trial = [];
event_marker = [];
conv_trial = DATA;

cond = event(:,8);
J_data = reshape(conv_trial,[size(conv_trial,1),size(conv_trial,2)*size(conv_trial,3)]);%将次个trial内所有神经元数据拼接为一行，每95bin为一个神经元随时间的spikecount情况
X = zscore(J_data); %数据标准化，效果好的情况与归一化相反
X = reshape(X,[size(conv_trial,1),size(conv_trial,2),size(conv_trial,3)]);
X = permute(X,[2 1 3]);X = reshape(X,[size(X,1),size(X,2)*size(X,3)]);
J_data = X;

clearvars -except data  J_data cond D conv_trial event_marker Nlist condition a con_t;
al_colour = [1 0 0; 1 0.5 0; 1 0 0.5; ...
    0.4 0.7 0.3; 0.2 0.8 0.1; 0.5 0.9 0.1; ...
    0 0.5 0.8; 0.58 0.8 1; 0.2 0.2 0.8; ...
    1 0.7 0; 1 0.9 0; 0.9 0.76 0.4; ...
    0 0.8 0.8; 0.2 0.8 0.8; 0.6 0.8 0.8; ...
    0.8 0 0.8; 0.8 0.2 0.8; 0.8 0.6 0.8];
%% draw CO figure
index_CO = cond<7;
clu_data_CO = J_data(index_CO,:);
clust_CO = cond(index_CO,:);
[coeff_CO,score_CO,~,~,explained_CO,mu] = pca(clu_data_CO);
cover_CO = cumsum(explained_CO);
n_PCs = find(cover_CO>90,1);%cover 90% 以上的pc数量

qerror_CO = [];
cverror_CO = [];
for i=1:n_PCs
    obj_CO = fitcdiscr(score_CO(:,1:i),clust_CO); %
    qerror_CO(i) = resubLoss(obj_CO);
    cvmodel = crossval(obj_CO,'kfold',10);
    cverror_CO(i) = kfoldLoss(cvmodel);
end
[~,n_PCs] = min(cverror_CO);
n_PCs = max(n_PCs,2);
obj_CO = fitcdiscr(score_CO(:,1:n_PCs),clust_CO); %
[w , ~] = eig(obj_CO.Sigma \ obj_CO.BetweenSigma);
[q , ~] = qr(w);

odr_CO = score_CO(:,1:n_PCs)*q(:,1:2);


colour = al_colour([1 4 7 10 13 16],:); %CO六方向
co = figure;
for i=1:length(unique(clust_CO))
    plot(odr_CO(clust_CO==i,1),odr_CO(clust_CO==i,2),'o','color',colour(i,:));
    p = ellipse(q(:,1:2)', cov(score_CO(clust_CO==i,1:n_PCs)), mean(score_CO(clust_CO==i,1:n_PCs))');
    line(p(1,:), p(2,:),'color',colour(i,:));
    hold on
end
[EVA] = cal_EVA(odr_CO,clu_data_CO);
legend(['Variance Explained ',num2str(EVA),'%']); legend('boxoff');
set(gca,'fontsize',10);
axis equal
set(gca,'ytick',[],'xtick',[],'xlim',1.3*xlim,'ylim',1.3*ylim);
box on

%  30 DR cond projected into CO space
score_all = bsxfun(@minus,J_data,mu)/coeff_CO'; %%reprojection test code MAY.6
odr_all = score_all(:,1:n_PCs)*q(:,1:2);
ss = rem([1:36],6);ss(ss==0)=6;
seq = [1:3:16 2:3:17 3:3:18];
colour = al_colour(seq,:); %in 18 condition order
mark  = {'o' 'd' '<' '*' '+' 's'};
figure;
for i=7:max(cond)
    plot(odr_all(cond==i,1),odr_all(cond==i,2),char(mark(floor(i/7)+1)),'color',colour(ss(i),:));
    p = ellipse(q(:,1:2)', cov(score_all(cond==i,1:n_PCs)), mean(score_all(cond==i,1:n_PCs))');
    line(p(1,:), p(2,:),'color',colour(ss(i),:));
    hold on
end
[EVA] = cal_EVA(odr_all(cond>6,:),J_data(cond>6,:));
legend(['Variance Explained ',num2str(EVA),'%']); legend('boxoff');
set(gca,'fontsize',10);
axis equal
set(gca,'ytick',[],'xtick',[],'xlim',1.3*xlim,'ylim',1.3*ylim);
box on

%% 6 directions with 6 types subplot

fd = rem(cond,6);fd(fd==0)=6;Clus_Dis_sou = [];
rng(2022);
for md = 1:6
    
    index = find(fd==md);
    clu_data = J_data(index,:);
    clust = cond(index,:);
    
    
    [~,score,~,~,explained,~] = pca(clu_data);
    [idx,~] = fscmrmr(score,clust);
    score = score(:,idx);
    cover = cumsum(explained);
    n_PCs = find(cover>90,1);
    
    % CROSS VALIDATION 寻找适合的参数数量
    qerror = [];
    cverror = [];
    for i=1:n_PCs
        obj = fitcdiscr(score(:,1:i),clust); %
        qerror(i) = resubLoss(obj);
        cvmodel = crossval(obj,'kfold',10);
        cverror(i) = kfoldLoss(cvmodel);
    end
    [M(md),n_PCs] = min(cverror);
    n_PCs = max(n_PCs,2);
    
    obj = fitcdiscr(score(:,1:n_PCs),clust);
    
    [w,~] = eig(obj.Sigma \ obj.BetweenSigma);
    q = gsog(w);
    
    odr = score(:,1:n_PCs)*q(:,1:2);
    
    %开始绘图

    mark  = {'o' 'd' '<' '*' '+' 's'};
    
    
    seq = unique(clust);
    [~] = cal_ma(odr,clust,seq,2);% 绘制马氏距离
    tmp = colour(md,:);
    mmp = [linspace(tmp(1),0.95,100)',linspace(tmp(2),0.95,100)',linspace(tmp(3),0.95,100)'];
    colormap(flipud(mmp))

    f=figure;
    for i= 1:6
        Clus_Dis_sou{md,i} = odr(clust==seq(i),:);
        a = scatter(odr(clust==seq(i),1),odr(clust==seq(i),2),char(mark(i)),'MarkerEdgeAlpha',.3,'MarkerEdgeColor',colour(ss(seq(i)),:));
        
        p = ellipse(q(:,1:2)', cov(score(clust==seq(i),1:n_PCs)), mean(score(clust==seq(i),1:n_PCs))');
        line(p(1,:), p(2,:),'color',colour(ss(seq(i)),:));
        hold on
        plot(mean(odr(clust==seq(i),1)),mean(odr(clust==seq(i),2)),char(mark(i)),'color',colour(ss(seq(i)),:),'MarkerSize',10.,'LineWidth',2);
    end
    axis equal
    set(gca,'ytick',[],'xtick',[],'xlim',1.3*xlim,'ylim',1.3*ylim);
    box off
    hold off
    ax=gca;ax.YAxis.Visible = 'off';   
    ax.XAxis.Visible = 'off';   
    set(gca,'color','none');
end
end