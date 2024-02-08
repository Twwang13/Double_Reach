% clear
% close all
% [filename, pathname] = uigetfile( {'*.*'}, 'Pick a file','on');
% [spike_processed,event] = process_data(pathname,filename);
function f = tt_analysis_plot(spike_processed,event,plot_name,zone) %,name
%%
% load(strcat(pathname,filename));
%
% condition = rec_beh(:,5);
% condition = condition(condition~=0);
%     condition(end) = [];
%
%
% for i = 1:length(trialTimeNum)-1
%     trial_rec_beh{i} = rec_beh(trialTimeNum(i):trialTimeNum(i+1)-1,:);
%     if any(trial_rec_beh{i}(:,3) == 1) && any(trial_rec_beh{i}(:,3) == 6)  && any(trial_rec_beh{i}(:,3) == 7) ; % && trial_rec_beh{i}(1,5)<7
%         correct_rec_beh{i} = trial_rec_beh{i};
%     else correct_rec_beh{i} = [];
%     end
% end
%
% condition(cellfun(@isempty,correct_rec_beh)) = [];
% correct_rec_beh(cellfun(@isempty,correct_rec_beh)) = [];
%
% move_onset = cellfun(@(x)x(x(:,3)==11 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% reach = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==12 & x(:,2)==0))-1)]*x(find(x(:,3)==12 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);%12有漏打现象
% ind2 = cellfun('length',reach)==0;
% ind = cellfun('length',move_onset)==0;%找出reach为空的trial，为CO，再从reach_two中对应的值转移到reach中，这样统一了CO的reach和DR的第一个reach
% ind3 = ind2-ind;
% correct_rec_beh = correct_rec_beh(ind3==0);
% condition = condition(ind3==0);%去掉reach漏打的trial
%
% % move_onset_two = cellfun(@(x)x(x(:,3)==13 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% % reach_two = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==14 & x(:,2)==0))-1)]*x(find(x(:,3)==14 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
% % ind2 = cellfun('length',reach_two)==0;
% % ind = cellfun('length',move_onset_two)==0;%找出reach为空的trial，为CO，再从reach_two中对应的值转移到reach中，这样统一了CO的reach和DR的第一个reach
% % ind3 = ind2-ind;
% % correct_rec_beh = correct_rec_beh(ind3==0);
% % condition = condition(ind3==0);%去掉reach_two漏打的trial
%
% move_onset = cellfun(@(x)x(x(:,3)==11 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% reach = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==12 & x(:,2)==0))-1)]*x(find(x(:,3)==12 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
% move_onset_two = cellfun(@(x)x(x(:,3)==13 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% reach_two = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==14 & x(:,2)==0))-1)]*x(find(x(:,3)==14 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
% touch = cellfun(@(x)x(x(:,3)==3 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% touch_two = cellfun(@(x)x(x(:,3)==7 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% ind = cellfun('length',move_onset)==0;%找出reach为空的trial，为CO，再从reach_two中对应的值转移到reach中，这样统一了CO的reach和DR的第一个reach
% move_onset(ind) = move_onset_two(ind);
% move_onset_two(ind) = [];
% reach(ind) = reach_two(ind);
% reach_two(ind) = [];
% touch(ind) = touch_two(ind);
% touch_two(ind) = [];
% que2 = cellfun('length',move_onset)==0;%从这里开始，是用于检查CO的是否漏打的代码，还需要重复跑上面的代码；
% que = cellfun('length',reach)==0;
% que = que + que2;%move 和 reach 仍有漏打
% ind4 = que==0;
% correct_rec_beh = correct_rec_beh(ind4);
% condition = condition(ind4);%去掉reach漏打的trial
% move_onset = cellfun(@(x)x(x(:,3)==11 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% reach = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==12 & x(:,2)==0))-1)]*x(find(x(:,3)==12 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
% move_onset_two = cellfun(@(x)x(x(:,3)==13 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% reach_two = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==14 & x(:,2)==0))-1)]*x(find(x(:,3)==14 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
% touch = cellfun(@(x)x(x(:,3)==3 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% touch_two = cellfun(@(x)x(x(:,3)==7 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
% ind = cellfun('length',move_onset)==0;%找出reach为空的trial，为CO，再从reach_two中对应的值转移到reach中，这样统一了CO的reach和DR的第一个reach
% move_onset(ind) = move_onset_two(ind);
% move_onset_two(ind) = [];
% reach(ind) = reach_two(ind);
% reach_two(ind) = [];
% touch(ind) = touch_two(ind);
% touch_two(ind) = [];
% % a=[1,zeros(1,length(find(correct_rec_beh{1,2}(:,3)==14 & correct_rec_beh{1,2}(:,2)==0))-1)]
% % b=correct_rec_beh{1,2}(find(correct_rec_beh{1,2}(:,3)==14 & correct_rec_beh{1,2}(:,2)==0),:)
% % a*b
% % c=cell2mat(c) 以上一句提取0 14 0,0 12 0marker唯一化的思路
% correct_rec_beh = cell2mat(correct_rec_beh');
%
% START = correct_rec_beh(correct_rec_beh(:,3)>100);
% END = correct_rec_beh(correct_rec_beh(:,3)==7);
% target_on = correct_rec_beh(correct_rec_beh(:,3)==1);
% go_signal = correct_rec_beh(correct_rec_beh(:,3)==6);
% move_onset = cell2mat(move_onset');
% reach = cell2mat(reach');
% touch = cell2mat(touch');
% move_onset_two = cell2mat(move_onset_two');
% reach_two = cell2mat(reach_two');
% touch_two = cell2mat(touch_two');
%
% event_trial_data = [START,target_on,go_signal,move_onset,reach,touch,END,condition];
% event_full = event_trial_data;
% event_full(event_trial_data(:,8)>6,9:11)=[move_onset_two,reach_two,touch_two];
% CO = event_full(event_trial_data(:,8)<7,:);
% DR = event_full(event_trial_data(:,8)>6,:);
% %CO = [START,target_on,go_signal,move_onset,reach,touch,END,condition]
% %DR = [START,target_on,go_signal,move_onset,reach,touch,END,condition,move_onset_two,reach_two,touch_two]
% %event_full = [START,target_on,go_signal,move_onset,reach,touch,END,condition,move_onset_two,reach_two,touch_two]
%
% clearvars -except CO DR event_trial_data spikeData pathname filename event_full;
%
% spikeTime = spikeData.times;
% binsize = 0.02;
% firstBin = [-0.6 -1 -0.4 -1];
% lastBin = [0.5 1 1 1];
% mBin = [0 0 -0.2 0 -0.08];
% nBin = [0.4 0.2 0 0.2 0];
%
% for i=1:length(spikeTime)%第一列是spike的时刻，第二列是对应的trial数，第三列是condition，第四列取bin marker信息
%     for k=1:length(event_trial_data)%数目都是一样的
%     if event_trial_data(k,1) < spikeTime(i,1) & event_trial_data(k,7) > spikeTime(i,1)
%         if event_trial_data(k,2)+mBin(1) < spikeTime(i,1) & event_trial_data(k,2)+nBin(1) > spikeTime(i,1)
%             spikeTime(i,4)=1;%第4列‘1’表示target_on
%         elseif event_trial_data(k,3)+mBin(2) < spikeTime(i,1) & event_trial_data(k,3)+nBin(2) > spikeTime(i,1)
%             spikeTime(i,5)=2;%第5列‘2’表示go_signal
%         elseif event_trial_data(k,4)+mBin(3) < spikeTime(i,1) & event_trial_data(k,4)+nBin(3) > spikeTime(i,1)
%             spikeTime(i,6)=3;%第6列‘3’表示move_onset
%         elseif event_trial_data(k,4)+mBin(4) < spikeTime(i,1) & event_trial_data(k,4)+nBin(4) > spikeTime(i,1)
%             spikeTime(i,7)=4;%第7列‘4’表示perimovement
%         elseif event_trial_data(k,5)+mBin(5) < spikeTime(i,1) & event_trial_data(k,5)+nBin(5) > spikeTime(i,1)
%             spikeTime(i,8)=5;%第8列‘5’表示reach
%         end
%        spikeTime(i,2)=k;
%        spikeTime(i,3)=event_trial_data(k,8);
%        trial_ind(k,:) = [k,event_trial_data(k,8)] ;
%
%     end
%     end
% end
% spikeTime(spikeTime(:,2)==0,:)=[];
% trial_ind(trial_ind(:,1)==0,:)=[];
% event = event_full(trial_ind(:,1),:);
% spike_processed = zeros(length(event),sum(spikeTime(:,2) == mode(spikeTime(:,2))));%找到单个trial中出现spike数最多的值，作为矩阵的2d
%
%
%
%     for k=1:length(trial_ind)
%         spike_processed(k,1:sum(spikeTime(:,2) == trial_ind(k))) = spikeTime(spikeTime(:,2) == trial_ind(k),1);
%     end
% spike_processed(spike_processed == 0) = nan; %至此数据整理完毕，主体为spike_processed,event,trial_ind,event_trial_data,DR
%%

binsize = 0.02;
firstBin = [-0.6 -1 -0.4 -1];
lastBin = [0.5 1 1 1];
mBin = [0 0 -0.2 0 -0.08];
nBin = [0.4 0.2 0 0.2 0];
ali_t = {'GO','MO','OFFSET'};

%下面开始画图
for p=3:5
    sequence = [9 3 2 6 12 13 ];
    spike_draw = spike_processed - event(:,p)*ones(1,size(spike_processed,2)); %对齐到哪里减去哪里
    
    f = figure (p) ;%(1)
    for i=1:6
        data = spike_draw(event(:,8)==i,:);
        %     c = mat2cell(data,ones(1,size(data,1)),[size(data,2)]);
        s(i) = subplot(3,5,sequence(i));
        %     [xPoints, yPoints] = SpikeRaster(c);
        %     psth(xPoints, yPoints, binsize, firstBin(2), lastBin(2),'k');
        psth_An(data,0.02,'k');
        xlim([-zone zone]);box off
        lim((i-1)*3+1) = max(get(gca,'Ylim'));
        hold on;%画CO
        
        data = spike_draw(event(:,8)==i+6,:);
        %     c = mat2cell(data,ones(1,size(data,1)),[size(data,2)]);
        s(i) = subplot(3,5,sequence(i));
        %     [xPoints, yPoints] = SpikeRaster(c);
        %     psth(xPoints, yPoints, binsize, firstBin(2), lastBin(2),'r');
        psth_An(data,0.02,'r');
        xlim([-zone zone]); box off
        lim((i-1)*3+2) = max(get(gca,'Ylim'));
        hold on;%画DR,CC
        
        data = spike_draw(event(:,8)==i+12,:);
        %     c = mat2cell(data,ones(1,size(data,1)),[size(data,2)]);
        s(i) = subplot(3,5,sequence(i));
        %     [xPoints, yPoints] = SpikeRaster(c);
        %     psth(xPoints, yPoints, binsize, firstBin(2), lastBin(2),'b');
        psth_An(data,0.02,'b');
        xlim([-zone zone]); box off
        lim(i*3) = max(get(gca,'Ylim'));
        hold on;%画DR,C
    end
    lim = max(lim);
    h = get(gcf,'Children');
    for i=1:6
        h(i).YLim = [0 lim];
    end
    for i=1:6
        data = spike_draw(event(:,8)==i,:);
        event_draw = event(event(:,8)==i,:);
        event_draw = event_draw(:,3:5) - event_draw(:,p)*ones(1,3);
        s(i) = subplot(3,5,sequence(i));
        rasters(data,'k',event_draw,p);
        hold on;%画CO
        
        data = spike_draw(event(:,8)==i+6,:);
        event_draw = event(event(:,8)==i+6,:);
        event_draw = event_draw(:,[3:5 9 10]) - event_draw(:,p)*ones(1,5);
        s(i) = subplot(3,5,sequence(i));
        rasters(data,'r',event_draw,p);
        hold on;%画DR,CC
        
        data = spike_draw(event(:,8)==i+12,:);
        event_draw = event(event(:,8)==i+12,:);
        event_draw = event_draw(:,[3:5 9 10]) - event_draw(:,p)*ones(1,5);
        s(i) = subplot(3,5,sequence(i));
        rasters(data,'b',event_draw,p);
        hold on;%画DR,C
    end
    
    draw_GO = spike_processed - event(:,3)*ones(1,size(spike_processed,2));
    draw_MO = spike_processed - event(:,4)*ones(1,size(spike_processed,2));
    draw_RE = spike_processed - event(:,5)*ones(1,size(spike_processed,2));
    peri={draw_GO,draw_MO,draw_RE};
    sequence_p = [5 10 15];
    clear mBin nBin;
    mBin = [-0.6 -0.2 -0.2];
    nBin = [0 0 0];
    for i = 1:3
        draw_data = peri{1,i};
        for k = 1:6
            data = draw_data(event(:,8)==k,:);
            c = mat2cell(data,ones(1,size(data,1)),[size(data,2)]);
            [xPoints, yPoints] = SpikeRaster(c);
            direction_ave_co(i,k) = direction_ave_cal(xPoints, yPoints, binsize, firstBin(2), lastBin(2), mBin(i), nBin(i));
            
            data = draw_data(event(:,8)==k+6,:);
            c = mat2cell(data,ones(1,size(data,1)),[size(data,2)]);
            [xPoints, yPoints] = SpikeRaster(c);
            direction_ave_ccw(i,k) = direction_ave_cal(xPoints, yPoints, binsize, firstBin(2), lastBin(2), mBin(i), nBin(i));
            
            data = draw_data(event(:,8)==k+12,:);
            c = mat2cell(data,ones(1,size(data,1)),[size(data,2)]);
            [xPoints, yPoints] = SpikeRaster(c);
            direction_ave_cw(i,k) = direction_ave_cal(xPoints, yPoints, binsize, firstBin(2), lastBin(2), mBin(i), nBin(i));
        end
        
        %     d0 = ceil(max(direction_ave_co)/5)*5;
        %     d1 = ceil(max(direction_ave_ccw)/5)*5;
        %     d2 = ceil(max(direction_ave_cw)/5)*5;
        %     r = max([d0 d1 d2]);
        
        %     s(i) = subplot(3,5,sequence_p(i));
        %     d0 = ceil(max(direction_ave_co)/5)*5;
        %     d1 = ceil(max(direction_ave_ccw)/5)*5;
        %     d2 = ceil(max(direction_ave_cw)/5)*5;
        %     r = max([d0 d1 d2]);
        %     if r == 0
        %         r = 5;
        %     end
        %     str = [num2str(r/2),' spike/s']
        %     plot([-1 -0.5],[-1 -1],'k','LineWidth',2);
        %     text(-1,-1.5,str,'FontSize',12);
        %     hold on
        %     pos = [-1 -1 2 2];
        %     rectangle('Position',pos,'Curvature',[1 1],'LineWidth',2);
        %     axis equal
        %     axis off
        %     hold on
        %     radius0 = direction_ave_co/r;
        %     radius1 = direction_ave_ccw/r;
        %     radius2 = direction_ave_cw/r;
        %     for j = 1:length(radius0)
        %         x0(j) = radius0(j)*cos(pi/3*(j-1));
        %         y0(j) = radius0(j)*sin(pi/3*(j-1));
        %     end
        %     plot([x0,x0(1)],[y0,y0(1)],'k','LineWidth',1);
        %     hold on;
        %     for j = 1:length(radius1)
        %         x1(j) = radius1(j)*cos(pi/3*(j-1));
        %         y1(j) = radius1(j)*sin(pi/3*(j-1));
        %     end
        %     plot([x1,x1(1)],[y1,y1(1)],'r','LineWidth',1);
        %     hold on;
        %     for j = 1:length(radius2)
        %         x2(j) = radius2(j)*cos(pi/3*(j-1));
        %         y2(j) = radius2(j)*sin(pi/3*(j-1));
        %     end
        %     plot([x2,x2(1)],[y2,y2(1)],'b','LineWidth',1);
        %     hold on;
        %     scatter(0,0,'filled','k');
        
    end
    d0 = ceil(max(max(direction_ave_co))/5)*5;
    d1 = ceil(max(max(direction_ave_ccw))/5)*5;
    d2 = ceil(max(max(direction_ave_cw))/5)*5;
    r = max([d0 d1 d2]);
    for i=1:3
        s(i) = subplot(3,5,sequence_p(i));
        str = [num2str(r/2),'sp/s'];
        plot([-1.1 -0.6],[-1 -1],'k','LineWidth',2);
        text(-2.2,-0.7,str,'FontSize',12);
        hold on
        pos = [-1 -1 2 2];
        rectangle('Position',pos,'Curvature',[1 1],'LineWidth',2);
        axis equal
        axis off
        hold on
        radius0 = direction_ave_co(i,:)/r;
        radius1 = direction_ave_ccw(i,:)/r;
        radius2 = direction_ave_cw(i,:)/r;
        [pd_co] = pd_calc(direction_ave_co(i,:),1);
        [pd_ccw] = pd_calc(direction_ave_ccw(i,:),1);
        [pd_cw] = pd_calc(direction_ave_cw(i,:),1);
        for j = 1:length(radius0)
            x0(j) = radius0(j)*cos(pi/3*(j-1));
            y0(j) = radius0(j)*sin(pi/3*(j-1));
        end
        plot([x0,x0(1)],[y0,y0(1)],'k','LineWidth',1);
        plot(1.15*cos(pd_co),1.15*sin(pd_co),'.','markersize',15,'MarkerEdgeColor','k');
        hold on;
        for j = 1:length(radius1)
            x1(j) = radius1(j)*cos(pi/3*(j-1));
            y1(j) = radius1(j)*sin(pi/3*(j-1));
        end
        plot([x1,x1(1)],[y1,y1(1)],'r','LineWidth',1);
        plot(1.15*cos(pd_ccw),1.15*sin(pd_ccw),'.','markersize',15,'MarkerEdgeColor','r');
        hold on;
        for j = 1:length(radius2)
            x2(j) = radius2(j)*cos(pi/3*(j-1));
            y2(j) = radius2(j)*sin(pi/3*(j-1));
        end
        plot([x2,x2(1)],[y2,y2(1)],'b','LineWidth',1);
        plot(1.15*cos(pd_cw),1.15*sin(pd_cw),'.','markersize',15,'MarkerEdgeColor','b');
        hold on;
        scatter(0,0,'filled','k');
    end
    
    set(gcf,'position',[350 350 1000 570]);
    
%     if ~isempty(plot_name)
%     %      print(f,'-dpng',['E:\Users\Leafstream\Documents\ME\DATA\B&C_dir\all_processed_data\psth\' name '_' ali_t{p-2} '.png']);
%     print(f,'-dpng',[plot_name,'_p',num2str(p) '.png']);
%     close all
%     else
%     end
    
end
end