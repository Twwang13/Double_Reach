function [co_spm,dr_spm,t,co_sp,dr_sp] = prep_softn(NKT,bin_w,stp,event,AE)


% 100ms size(NKT,1) 2301ms MO-1.3-1
if isempty(bin_w) && isempty(stp)
    data = softnorm(NKT*1000);
    t=-1300:1000;
else
    bin = [];
    
    % aligned time to GO
    if  ~isempty(AE)
        GM = fix(1000*(event(:,4)-event(:,3)));
        Atime = [-1300:1000]+GM;
        X = [];
        for i = 1:size(NKT,2)
            p_point = [find(Atime(i,:)==0)-599,find(Atime(i,:)==0)];
            X(:,i,:) = NKT(:,i,p_point(1):p_point(2));
        end
        for i =1:fix((size(X,3)-bin_w)/stp+1)
            bin(:,:,i)=sum(X(:,:,stp*(i-1)+1:stp*(i-1)+bin_w),3);
        end
        t = (-600:stp:0-bin_w)+(bin_w/2);
    else
        
        for i =1:fix((length(NKT)-bin_w)/stp+1)
            bin(:,:,i)=sum(NKT(:,:,stp*(i-1)+1:stp*(i-1)+bin_w),3);
        end
        t = (-1300:stp:1000-bin_w)+(bin_w/2);
    end
    data=bin*1000/bin_w;   
    data = softnorm(data);
end

figure
plot(squeeze(mean(data,2))');

cond=event(:,8);
fdir = rem(cond,6);
fdir(fdir==0)=6;
f_ang = 60*(fdir-1);

co_a = f_ang(cond<7);
co_sp = data(:,cond<7,:);


%% 整理拟合数据
co_spm = zeros(size(data,1),6,size(co_sp,3));
for c = 1:6
    co_spm(:,c,:) = mean(data(:,cond==c,:),2);
    std_co(:,c,:) = std(data(:,cond==c,:),0,2);
    SE_co(:,c,:) = std_co(:,c,:)/sqrt(sum(cond==c));
end
co_am = 0:60:300;

dr_fa = f_ang(cond>6);     %第一运动角度
dr_sp = data(:,cond>6,:);
dr_con = cond(cond>6);
dr_ind(dr_con<13)=-1; %120cccw
dr_ind((dr_con>12) & (dr_con<19))=1; %120cw
dr_ind((dr_con>18) & (dr_con<25))=0; %180
dr_ind((dr_con>24) & (dr_con<31))=-2; %60ccw
dr_ind(dr_con>30) = 2;
dr_sa = dr_fa+60*dr_ind'-180;%中心坐标
dr_sma = dr_fa+30*dr_ind'-180;%未来手坐标
dr_spm = zeros(size(data,1),12,size(dr_sp,3));dr_sap=[];

for c = 1:30
    dr_spm(:,c,:) = mean(dr_sp(:,dr_con==c+6,:),2);
    %     var_dr(:,c,:) = var(dr_sp(:,dr_con==c+6,:),1,2);
    std_dr(:,c,:) = std(dr_sp(:,dr_con==c+6,:),0,2);
    SE_dr(:,c,:) = std_dr(:,c,:)/sqrt(sum(dr_con==c+6));
    dr_am(c) = mean(dr_fa(dr_con==c+6)); %第一伸手方向
    dr_asm(c) = mean(dr_sma(dr_con==c+6));%手坐标第二点
    dr_sam(c) = mean(dr_sa(dr_con==c+6));%中心坐标第二点
    dr_sap(:,c,:) = mean(dr_sp(:,dr_sma==dr_asm(c),:),2);
end
% soft normalization
% tmp = reshape(co_spm,[size(co_spm,1),size(co_spm,2)*size(co_spm,3)]);
% SE_co_N = SE_co./repmat((range(tmp')'+5),[1 size(SE_co,2) size(SE_co,3)]);
% tmp = tmp./repmat(range(tmp')'+5,[1 size(tmp,2)]);
% xx = reshape(tmp,[size(co_spm,1),size(co_spm,2),size(co_spm,3)]);
% co_spm = xx;
% co_sp = softnorm(co_sp);
%
%
% tmp = reshape(dr_spm,[size(dr_spm,1),size(dr_spm,2)*size(dr_spm,3)]);
% % var_dr_N = var_dr./repmat((range(tmp')'+5).^2,[1 size(var_dr,2) size(var_dr,3)]);
% % std_dr_N = std_dr./repmat((range(tmp')'+5),[1 size(std_dr,2) size(std_dr,3)]);
% SE_dr_N = SE_dr./repmat((range(tmp')'+5),[1 size(SE_dr,2) size(SE_dr,3)]);
% tmp = tmp./repmat(range(tmp')'+5,[1 size(tmp,2)]); %softmax
% xx = reshape(tmp,[size(dr_spm,1),size(dr_spm,2),size(dr_spm,3)]);
% dr_spm = xx;
% dr_sp = softnorm(dr_sp);
end




