function corval = plotHandtrj_emg
load('data/behData.mat');
%% cut into single trials
tri_data=[];
for tri = 1:length(event_full)
    if event_full(tri,10)==0
        tri_data{tri} = use_data(vicon_time>(event_full(tri,2)) & vicon_time<(event_full(tri,5)),:);
    else
        tri_data{tri} = use_data(vicon_time>(event_full(tri,2)) & vicon_time<(event_full(tri,10)),:);
    end
end

%% resample
kao = [];
for i = 1:length(tri_data)
    play = tri_data{i};
    x = 1:length(play);
    t = linspace(1, length(play), 200);
    kao(:,:,i) = pchip(x,play',t)';
end
tri_data = kao;

%% plot Fig 1b
al_colour = [1 0 0; 1 0.5 0; 1 0 0.5; ...
    0.4 0.7 0.3; 0.2 0.8 0.1; 0.5 0.9 0.1; ...
    0 0.5 0.8; 0.58 0.8 1; 0.2 0.2 0.8; ...
    1 0.7 0; 1 0.9 0; 0.9 0.76 0.4; ...
    0 0.8 0.8; 0.2 0.8 0.8; 0.6 0.8 0.8; ...
    0.8 0 0.8; 0.8 0.2 0.8; 0.8 0.6 0.8];
seq = [1:3:16 2:3:17 3:3:18];
colour = al_colour(seq,:); %in 18 condition order

figure
for d = 1:6
    for i = 0:2
        plot_Data = tri_data(:,:,event_full(:,8)==d+i*6);
        plot_Data = plot_Data(:,:,2:end);
        subplot(3,2,d)
        for tri = 1:size(plot_Data,3)
            %             plot3(plot_Data(:,1,tri),plot_Data(:,2,tri),plot_Data(:,3,tri),'color',colour(d+i*6,:));
            hold on
            %             plot3(plot_Data(tri).hand(1,1),plot_Data(tri).hand(1,2),plot_Data(tri).hand(1,3),'color',colour(d+i*6,:),'marker','o');
            plot(-plot_Data(:,1,tri),plot_Data(:,3,tri),'color',colour(d+i*6,:));
            hold on
            %             plot(-plot_Data(tri).hand(1,1),plot_Data(tri).hand(1,3),'color',colour(d+i*6,:),'marker','o');
            box off
            xlim([-300 300]);ylim([-50 550]); axis square;ax=gca;ax.YAxis.Visible = 'off';   % 设置y轴不可见
            ax.XAxis.Visible = 'off';   % 设置y轴不可见
        end
    end
end
set(gcf,'position',[1000 380 300 445]);

%% plot Fig1c
condition = event_full(:,8);

%% speed profile
% 数据采样按event前后固定时间长度采样
ali = 4; %align to screen MO
t_range = [0.475 0.775];
tri_data = zeros(sum(t_range)*100,size(use_data,2),length(event_full));
for tri = 1:length(event_full)
    tri_data(:,:,tri) = use_data(vicon_time>(event_full(tri,ali)-t_range(1)) & vicon_time<(event_full(tri,ali)+t_range(2)),:);
end

speed = [];
for tri = 1:length(tri_data)
    plot_Data = tri_data(:,:,tri);
    speed(tri,:) = sum(abs(plot_Data(:,4:6)).^2,2).^(1/2)';
    %     plot(speed,'color',colour(event_full(tri,8),:));
    %     hold on
end

% norm + corr acc get speed and condtion and event_full
% 10ms bin 125bins from 475-MO-775
figure
t_end=58; % avg offset
t = -0.474:0.01:0.775;
sequence = [9 3 2 6 12 13 ];
for d  =1:6
    
    s=subplot(3,5,sequence(d));
    vel_sr=nanmean(speed(condition==d,:));
    vel_ccw=nanmean(speed(condition==d+6,:));
    vel_cw=nanmean(speed(condition==d+12,:));
    plot(t,vel_sr,'k');hold on
    E=std(speed(condition==d,:))/sqrt(sum(condition==d));
    fill([t fliplr(t)],[vel_sr+2*E fliplr(vel_sr-2*E)],'k','LineStyle','none');alpha(.3);

    plot(t,vel_ccw,'r');hold on
    E=std(speed(condition==d+6,:))/sqrt(sum(condition==d+6));
    fill([t fliplr(t)],[vel_ccw+2*E fliplr(vel_ccw-2*E)],'r','LineStyle','none');alpha(.3);

    plot(t,vel_cw,'b');hold on
    E=std(speed(condition==d+12,:))/sqrt(sum(condition==d+12));
    fill([t fliplr(t)],[vel_cw+2*E fliplr(vel_cw-2*E)],'b','LineStyle','none');alpha(.3);

    ylim([-2000,4500])
    plot(-0.15,-100,'k.','markersize',12);
    a = ylim;
    plot(t(t_end),-100,'k.','markersize',12);
    xlim([-0.475 0.775])
    cor_vel(d) = corr(vel_sr(1:t_end)',mean([vel_ccw(1:t_end);vel_cw(1:t_end)])');

    box off
    YAxis.Visible = 'off';XAxis.Visible = 'off';
    set(gca,'Visible','off');
    set(gcf,'position',[679 127 613 420]);
    ylim([-2000,4500])
end
corval = mean(cor_vel);
std(cor_vel);

%% plot emg
sensornum = fix(size(emg,2)/20);
lab_s = event_full(condition<7,:);
lab_d = event_full(condition>6,:);

for sensor = 1:sensornum
    emg_t=emg(:,(sensor-1)*20+1);
    emg_v=emg(:,(sensor-1)*20+2); %in volt，每个传感器数据为10D，EMG+XYZ+pitch,yaw,roll+magXYZ,每个参量前跟随对应时间，故为20D
    acct=emg(:,(sensor-1)*20+3);          %acc 时间，所有惯性传感器共享
    accxyz=emg(:,[(sensor-1)*20+4 (sensor-1)*20+6 (sensor-1)*20+8]); %in xyz order 当sensor1固定在末端时，用468来描绘endpoint的运动参数
    mag = emg(:,[(sensor-1)*20+10 (sensor-1)*20+12 (sensor-1)*20+14]);
    magt = emg(:,(sensor-1)*20+9);
    
    %processing single reach
    for i=1:length(lab_s)
        emg_s_ind{i} = find(emg_t>lab_s(i,1) & emg_t<lab_s(i,7));
        emg_s_t{i} = emg_t(emg_s_ind{i})-lab_s(i,4); %减去lab 2就是对齐到2列的event
        emg_s_ct{i} = emg_t(emg_s_ind{i})-lab_s(i,3);%对齐到go开始
        emg_s_v{i} = emg_v(emg_s_ind{i});
        accxyz_s{i} = accxyz(acct>lab_s(i,1) & acct<lab_s(i,7),:);
        acct_s{i} = acct(acct>lab_s(i,1) & acct<lab_s(i,7))-lab_s(i,4);        % time aligned to MO
        acct_sc{i} = acct(acct>lab_s(i,1) & acct<lab_s(i,7))-lab_s(i,3);%*
        mag_ts{i} = magt(magt>lab_s(i,1) & magt<lab_s(i,7))-lab_s(i,3);
        magxyz{i}=mag(magt>lab_s(i,1) & magt<lab_s(i,7),:);
    end
    %processing double reach
    for i=1:length(lab_d)
        emg_d_ind{i} = find(emg_t>lab_d(i,1) & emg_t<lab_d(i,7));
        emg_d_t{i} = emg_t(emg_d_ind{i})-lab_d(i,4);
        emg_d_ct{i} = emg_t(emg_d_ind{i})-lab_d(i,3);%*
        emg_d_v{i} = emg_v(emg_d_ind{i});
        accxyz_d{i} = accxyz(acct>lab_d(i,1) & acct<lab_d(i,7),:);
        acct_d{i} = acct(acct>lab_d(i,1) & acct<lab_d(i,7))-lab_d(i,4);        % time aligned to MO
        acct_dc{i} = acct(acct>lab_d(i,1) & acct<lab_d(i,7))-lab_d(i,3);
    end
    
    %produce single reach matrix
    for j = 1:length(emg_s_t)
        ind_s_t = emg_s_t{j};
        ind_s_ct = emg_s_ct{j}; %*
        ind_s_v = emg_s_v{j};
        ind_s_a = accxyz_s{j};
        ind_s_at = acct_s{j};
        ind_sc_at = acct_sc{j}; %*
        m = 1;
        n=1;
        for i = -0.475:0.025:0.775 %bin is 50ms, step is 25ms
            bin_t(m) = i;
            bin_s_v(j,m) = rms(ind_s_v(ind_s_t > i-0.025 & ind_s_t < i+0.025));
            %     bin_s_a(j,m) = mean(ind_s_a(ind_s_at > i-0.025 & ind_s_at < i+0.025));
            bin_s_axyz(j,m,1:3) = nanmean(ind_s_a(ind_s_at > i-0.025 & ind_s_at < i+0.025,:)); %更改
            m=m+1;
        end
        for i = -0.875:0.025:0.175 %*
            bin_c_t(n) = i;
            bin_s_vc(j,n) = rms(ind_s_v(ind_s_ct > i-0.025 & ind_s_ct < i+0.025));
            bin_sc_axyz(j,n,1:3) = nanmean(ind_s_a(ind_sc_at > i-0.025 & ind_sc_at < i+0.025,:));
            n=n+1;
        end
    end
    bin_all_t = [bin_c_t,nan,bin_t+0.7];
    bin_s_v_all = [bin_s_vc,zeros(size(bin_s_v,1),1)*nan,bin_s_v];
    bin_s_axyz_all = [bin_sc_axyz,zeros(size(bin_s_axyz,1),1,size(bin_s_axyz,3))*nan,bin_s_axyz];
    
    %produce double reach matrix
    for j = 1:length(emg_d_t)
        ind_d_t = emg_d_t{j};
        ind_d_ct = emg_d_ct{j};
        ind_d_v = emg_d_v{j};
        ind_d_a = accxyz_d{j};
        ind_d_at = acct_d{j};
        ind_dc_at = acct_dc{j};
        m = 1;
        n = 1;
        for i = -0.475:0.025:0.775 %bin is 50ms, step is 25ms
            bin_d_v(j,m) = rms(ind_d_v(ind_d_t > i-0.025 & ind_d_t < i+0.025));
            %     bin_d_a(j,m) = mean(ind_d_a(ind_d_at > i-0.025 & ind_d_at < i+0.025));
            bin_d_axyz(j,m,1:3) = nanmean(ind_d_a(ind_d_at > i-0.025 & ind_d_at < i+0.025,:));
            m=m+1;
        end
        for i = -0.875:0.025:0.175 %*
            bin_d_vc(j,n) = rms(ind_d_v(ind_d_ct > i-0.025 & ind_d_ct < i+0.025));
            bin_dc_axyz(j,n,1:3) = nanmean(ind_d_a(ind_dc_at > i-0.025 & ind_dc_at < i+0.025,:));
            n=n+1;
        end
    end
    bin_d_v_all = [bin_d_vc,zeros(size(bin_d_v,1),1)*nan,bin_d_v];
    bin_d_axyz_all = [bin_dc_axyz,zeros(size(bin_d_axyz,1),1,size(bin_d_axyz,3))*nan,bin_d_axyz];
    
    %%%  new 绘图
    L_emg = cat(1,bin_s_v_all,bin_d_v_all);
    labels = cat(1,lab_s(:,8),lab_d(:,8));
    t_end = 68;
    up = 0.0035;
    down = min(min(L_emg));
    figure
    for d  =1:6
        s=subplot(3,5,sequence(d));
        sr=mean(L_emg(labels==d,:));ccw=mean(L_emg(labels==d+6,:));cw=mean(L_emg(labels==d+12,:));
        plot(bin_all_t,sr,'k');hold on
        E=std(L_emg(labels==d,:))/sqrt(sum(labels==d));
        fill([bin_all_t(1:length(bin_c_t)) fliplr(bin_all_t(1:length(bin_c_t)))],[sr(1:length(bin_c_t))+2*E(1:length(bin_c_t)) fliplr(sr(1:length(bin_c_t))-2*E(1:length(bin_c_t)))],'k','LineStyle','none');
        alpha(.3);
        fill([bin_all_t(length(bin_c_t)+2:end) fliplr(bin_all_t(length(bin_c_t)+2:end))],[sr(length(bin_c_t)+2:end)+2*E(length(bin_c_t)+2:end) fliplr(sr(length(bin_c_t)+2:end)-2*E(length(bin_c_t)+2:end))],'k','LineStyle','none');
        alpha(.3);
        %         fill([bin_all_t fliplr(bin_all_t)],[sr+2*E fliplr(sr-2*E)],'k','LineStyle','none');alpha(.3);
        
        plot(bin_all_t,ccw,'r');hold on
        E=std(L_emg(labels==d+6,:))/sqrt(sum(labels==d+6));
        fill([bin_all_t(1:length(bin_c_t)) fliplr(bin_all_t(1:length(bin_c_t)))],[ccw(1:length(bin_c_t))+2*E(1:length(bin_c_t)) fliplr(ccw(1:length(bin_c_t))-2*E(1:length(bin_c_t)))],'r','LineStyle','none');
        alpha(.3);
        fill([bin_all_t(length(bin_c_t)+2:end) fliplr(bin_all_t(length(bin_c_t)+2:end))],[ccw(length(bin_c_t)+2:end)+2*E(length(bin_c_t)+2:end) fliplr(ccw(length(bin_c_t)+2:end)-2*E(length(bin_c_t)+2:end))],'r','LineStyle','none');
        alpha(.3);
     
        plot(bin_all_t,cw,'b');hold on
        E=std(L_emg(labels==d+12,:))/sqrt(sum(labels==d+12));
        fill([bin_all_t(1:length(bin_c_t)) fliplr(bin_all_t(1:length(bin_c_t)))],[cw(1:length(bin_c_t))+2*E(1:length(bin_c_t)) fliplr(cw(1:length(bin_c_t))-2*E(1:length(bin_c_t)))],'b','LineStyle','none');
        alpha(.3);
        fill([bin_all_t(length(bin_c_t)+2:end) fliplr(bin_all_t(length(bin_c_t)+2:end))],[cw(length(bin_c_t)+2:end)+2*E(length(bin_c_t)+2:end) fliplr(cw(length(bin_c_t)+2:end)-2*E(length(bin_c_t)+2:end))],'b','LineStyle','none');
        alpha(.3);
 
        a = [down-up/10 up];
        ylim(a)
        plot(0.7-0.15,a(1),'k.','markersize',12);%补偿屏幕延时的MO
        plot(bin_all_t(t_end),a(1),'k.','markersize',12);%ME
        plot(0,a(1),'k.','markersize',12);%GO
        xlim([-0.85 1.475]);
        box off
        YAxis.Visible = 'off';XAxis.Visible = 'off';
        set(gca,'Visible','off');
        set(gcf,'position',[679 127 613 420]);


    end
    plot([1.275 1.475],[a(1) a(1)],'color',[150 150 150]/255,'linewidth',3);
end
