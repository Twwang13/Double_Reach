%% preprocessing
load('dataNKT_158.mat');

% soft-normalized each neuron FR
[co_spm,dr_spm,t,co_sp,dr_sp] = prep_softn(NKT(41:end,:,:),[],[],event,[]); % n*c*t avgfr

% mean-centered firing rate within each set of conditions
d_co = co_spm-mean([co_spm],2);d_cosp = co_sp-mean([co_spm,dr_spm],2);
d_dr = dr_spm-mean([dr_spm],2);d_drsp = dr_sp-mean([co_spm,dr_spm],2);
d_all = [d_co,d_dr];data = [d_cosp,d_drsp];
cond = event(:,8);
ind = [find(cond<7);find(cond>6)];cond = cond(ind);

J_data=data;
% J_data cond
J_data = permute(J_data,[2,1,3]);

fdir = rem(cond,6);
fdir(fdir==0)=6;
f_ang = 60*(fdir-1);  % 第一运动角度

dr_fa = f_ang(cond>6);     %第一运动角度
dr_data = J_data(cond>6,:,:);
dr_con = cond(cond>6);
dr_ind(dr_con<13)=-1; %120cccw
dr_ind((dr_con>12) & (dr_con<19))=1; %120cw
dr_ind((dr_con>18) & (dr_con<25))=0; %180
dr_ind((dr_con>24) & (dr_con<31))=-2; %60ccw
dr_ind(dr_con>30) = 2;
dr_sa = dr_fa+60*dr_ind'-180;%中心坐标
dr_sma = dr_fa+30*dr_ind'-180;%未来手坐标

% dr_sa dr_sma  dr_data dr_fa
% 取 bin
bin_w = 300; 
stp = 20;
bin=[];
for i =1:fix((size(dr_data,3)-bin_w)/stp+1)
    bin(:,:,i)=sum(dr_data(:,:,stp*(i-1)+1:stp*(i-1)+bin_w),3);
end
t = (-1300:stp:1000-bin_w)+(bin_w/2)+150;

% 解码
cverror_f=[];cverror_fp=[];cverror_sma=[];cverror_sp=[];
for T = 1:100

for i = 1:size(bin,3)
    obj = fitcdiscr(bin(:,:,i),dr_fa); %fa
    cvmodel = crossval(obj,'kfold',10);
    cverror_f(T,i) = kfoldLoss(cvmodel);
    
    obj = fitcdiscr(bin(:,:,i),dr_fa(randperm(length(dr_fa)))); %fa shuffle
    cvmodel = crossval(obj,'kfold',10);
    cverror_fp(T,i) = kfoldLoss(cvmodel);
    
    obj = fitcdiscr(bin(:,:,i),dr_sma); %sma
    cvmodel = crossval(obj,'kfold',10);
    cverror_sma(T,i) = kfoldLoss(cvmodel);
    
    obj = fitcdiscr(bin(:,:,i),dr_sma(randperm(length(dr_sma)))); %sma
    cvmodel = crossval(obj,'kfold',10);
    cverror_sp(T,i) = kfoldLoss(cvmodel);
    
%     obj = fitcdiscr(bin(:,:,i),dr_sa); %sa
%     cvmodel = crossval(obj,'kfold',10);
%     cverror_sa(i) = kfoldLoss(cvmodel);
end
end
figure
plot(t,1-mean(cverror_f),'k','LineWidth',1.5);
hold on
h1 = fill([t fliplr(t)],[1-mean(cverror_f)+std(1-cverror_f) fliplr(1-mean(cverror_f)-std(1-cverror_f))],'k');
alpha(.3);
set(h1,'LineStyle','none');
plot(t,1-mean(cverror_sma),'r','LineWidth',1.5);
h2 = fill([t fliplr(t)],[1-mean(cverror_sma)+std(1-cverror_sma) fliplr(1-mean(cverror_sma)-std(1-cverror_sma))],'r');
alpha(.3);
set(h2,'LineStyle','none');
% plot(t,1-cverror_sa,'b');
% plot(t,1-mean(cverror_sp),'g','LineWidth',1.5);
h3 = fill([t fliplr(t)],[1-mean(cverror_sp)+std(1-cverror_sp) fliplr(1-mean(cverror_sp)-std(1-cverror_sp))],'r');
alpha(.2);
set(h3,'LineStyle','none');
h4 = fill([t fliplr(t)],[1-mean(cverror_fp)+std(1-cverror_fp) fliplr(1-mean(cverror_fp)-std(1-cverror_fp))],'k');
alpha(.2);
set(h4,'LineStyle','none');

xlim([t(1),t(end)]);
ylim([0,1])
box off




















