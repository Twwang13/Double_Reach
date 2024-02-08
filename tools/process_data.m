function [spike_processed,event] = process_data(pathname,filename)

load(strcat(pathname,filename));

%处理spta用，对齐肌电数据

spikeData.times = spikeData.times(:)-rec_beh(1,1);
rec_beh(:,1) = rec_beh(:,1)-rec_beh(1,1);

condition = rec_beh(:,5);
condition = condition(condition~=0);
condition(end) = [];


for i = 1:length(trialTimeNum)-1
    trial_rec_beh{i} = rec_beh(trialTimeNum(i):trialTimeNum(i+1)-1,:);
    if any(trial_rec_beh{i}(:,3) == 1) && any(trial_rec_beh{i}(:,3) == 6)  && any(trial_rec_beh{i}(:,3) == 7) ; % && trial_rec_beh{i}(1,5)<7
        correct_rec_beh{i} = trial_rec_beh{i};
    else correct_rec_beh{i} = [];
    end
end

condition(cellfun(@isempty,correct_rec_beh)) = [];
correct_rec_beh(cellfun(@isempty,correct_rec_beh)) = [];

move_onset = cellfun(@(x)x(x(:,3)==11 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
reach = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==12 & x(:,2)==0))-1)]*x(find(x(:,3)==12 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);%12有漏打现象
ind2 = cellfun('length',reach)==0;
ind = cellfun('length',move_onset)==0;%找出reach为空的trial，为CO，再从reach_two中对应的值转移到reach中，这样统一了CO的reach和DR的第一个reach
ind3 = ind2-ind;
correct_rec_beh = correct_rec_beh(ind3==0);
condition = condition(ind3==0);%去掉reach漏打的trial

move_onset_two = cellfun(@(x)x(x(:,3)==13 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
reach_two = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==14 & x(:,2)==0))-1)]*x(find(x(:,3)==14 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
ind2 = cellfun('length',reach_two)==0;
ind = cellfun('length',move_onset_two)==0;%找出reach为空的trial，为CO，再从reach_two中对应的值转移到reach中，这样统一了CO的reach和DR的第一个reach
ind3 = ind2-ind;
correct_rec_beh = correct_rec_beh(ind3==0);
condition = condition(ind3==0);%去掉reach_two漏打的trial

move_onset = cellfun(@(x)x(x(:,3)==11 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
move_onset_two = cellfun(@(x)x(x(:,3)==13 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
reach = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==12 & x(:,2)==0))-1)]*x(find(x(:,3)==12 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
reach_two = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==14 & x(:,2)==0))-1)]*x(find(x(:,3)==14 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
touch = cellfun(@(x)x(x(:,3)==3 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
touch_two = cellfun(@(x)x(x(:,3)==7 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
ind = cellfun('length',reach)==0;%找出reach为空的trial，为CO，再从reach_two中对应的值转移到reach中，这样统一了CO的reach和DR的第一个reach
move_onset(ind) = move_onset_two(ind);
move_onset_two(ind) = [];
reach(ind) = reach_two(ind);
reach_two(ind) = [];
touch(ind) = touch_two(ind);
touch_two(ind) = [];
% a=[1,zeros(1,length(find(correct_rec_beh{1,2}(:,3)==14 & correct_rec_beh{1,2}(:,2)==0))-1)]
% b=correct_rec_beh{1,2}(find(correct_rec_beh{1,2}(:,3)==14 & correct_rec_beh{1,2}(:,2)==0),:)
% a*b
% c=cell2mat(c) 以上一句提取0 14 0,0 12 0marker唯一化的思路
correct_rec_beh = cell2mat(correct_rec_beh');

START = correct_rec_beh(correct_rec_beh(:,3)>100);
END = correct_rec_beh(correct_rec_beh(:,3)==7);
target_on = correct_rec_beh(correct_rec_beh(:,3)==1);
go_signal = correct_rec_beh(correct_rec_beh(:,3)==6);
move_onset = cell2mat(move_onset');
reach = cell2mat(reach');
touch = cell2mat(touch');
move_onset_two = cell2mat(move_onset_two');
reach_two = cell2mat(reach_two');
touch_two = cell2mat(touch_two');

event_trial_data = [START,target_on,go_signal,move_onset,reach,touch,END,condition];
event_full = event_trial_data;
event_full(event_trial_data(:,8)>6,9:11)=[move_onset_two,reach_two,touch_two];
CO = event_full(event_trial_data(:,8)<7,:);
DR = event_full(event_trial_data(:,8)>6,:);
%CO = [START,target_on,go_signal,move_onset,reach,touch,END,condition]
%DR = [START,target_on,go_signal,move_onset,reach,touch,END,condition,move_onset_two,reach_two,touch_two]
%event_full = [START,target_on,go_signal,move_onset,reach,touch,END,condition,move_onset_two,reach_two,touch_two]

clearvars -except CO DR event_trial_data spikeData pathname filename event_full;

spikeTime = spikeData.times;
% binsize = 0.02;
% firstBin = [-0.6 -1 -0.4 -1]; 
% lastBin = [0.5 1 1 1];
mBin = [0 0 -0.2 0 -0.08];
nBin = [0.4 0.2 0 0.2 0];

for i=1:length(spikeTime)%第一列是spike的时刻，第二列是对应的trial数，第三列是condition，第四列取bin marker信息
    for k=1:length(event_trial_data)%数目都是一样的
    if event_trial_data(k,1) < spikeTime(i,1) & event_trial_data(k,7) > spikeTime(i,1)
        if event_trial_data(k,2)+mBin(1) < spikeTime(i,1) & event_trial_data(k,2)+nBin(1) > spikeTime(i,1)
            spikeTime(i,4)=1;%第4列‘1’表示target_on
        elseif event_trial_data(k,3)+mBin(2) < spikeTime(i,1) & event_trial_data(k,3)+nBin(2) > spikeTime(i,1)
            spikeTime(i,5)=2;%第5列‘2’表示go_signal
        elseif event_trial_data(k,4)+mBin(3) < spikeTime(i,1) & event_trial_data(k,4)+nBin(3) > spikeTime(i,1)
            spikeTime(i,6)=3;%第6列‘3’表示move_onset 
        elseif event_trial_data(k,4)+mBin(4) < spikeTime(i,1) & event_trial_data(k,4)+nBin(4) > spikeTime(i,1)
            spikeTime(i,7)=4;%第7列‘4’表示perimovement
        elseif event_trial_data(k,5)+mBin(5) < spikeTime(i,1) & event_trial_data(k,5)+nBin(5) > spikeTime(i,1)
            spikeTime(i,8)=5;%第8列‘5’表示reach
        end            
       spikeTime(i,2)=k;
       spikeTime(i,3)=event_trial_data(k,8);
       trial_ind(k,:) = [k,event_trial_data(k,8)] ;
    
    end
    end
end

spikeTime(spikeTime(:,2)==0,:)=[];
trial_ind(trial_ind(:,1)==0,:)=[];
event = event_full(trial_ind(:,1),:);
spike_processed = zeros(size(event,1),sum(spikeTime(:,2) == mode(spikeTime(:,2))));



    for k=1:length(event)
        spike_processed(k,1:sum(spikeTime(:,2) == trial_ind(k,1))) = spikeTime(spikeTime(:,2) == trial_ind(k,1),1);
    end
spike_processed(spike_processed == 0) = nan; %至此数据整理完毕，主体为spike_processed,event,trial_ind,event_trial_data,DR
end

