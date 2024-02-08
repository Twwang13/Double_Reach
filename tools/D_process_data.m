function [spike_processed,event] = D_process_data(pathname,filename)

load(strcat(pathname,filename));
start_index = find(rec_beh(:,2)==9);
END_index = find(rec_beh(:,2)==18);
sum(END_index)-sum(start_index-1)==END_index(end); %检验是否匹配和多打
start = rec_beh(start_index,:);
END = rec_beh(END_index,:);
ee = who('-regexp','error');
eval(['errorData=', ee{1},';']);

trial_rec_beh = [];
correct_rec_beh = [];
for i=1:length(END)
    trial_rec_beh{i} = rec_beh(start_index(i):END_index(i),:);
    if errorData(i) == 0;
        correct_rec_beh{i} = trial_rec_beh{i};
    else correct_rec_beh{i} = [];
    end
end
condition = condition(errorData==0);
correct_rec_beh(cellfun(@isempty,correct_rec_beh)) = [];
event_rec_data = cell2mat(correct_rec_beh');
cue_on = event_rec_data(event_rec_data(:,2)==3,1);
go_signal = event_rec_data(event_rec_data(:,2)==5,1);
MO = event_rec_data(event_rec_data(:,2)==6,1);
reach = event_rec_data(event_rec_data(:,2)==14 | event_rec_data(:,2)==7,1);
touch = event_rec_data(event_rec_data(:,2)==15 | event_rec_data(:,2)==8,1);
event_full = zeros(length(condition),11);
event_full(:,1) = start(errorData==0);
event_full(:,2) = cue_on;
event_full(:,3) = go_signal;
event_full(:,4) = MO;
event_full(:,5) = reach;
event_full(:,6) = touch;
event_full(:,7) = END(errorData==0);
event_full(:,8) = condition;
event_full(condition>6,9) = event_rec_data(event_rec_data(:,2)==10,1);
event_full(condition>6,10) = event_rec_data(event_rec_data(:,2)==11,1);
event_full(condition>6,11) = event_rec_data(event_rec_data(:,2)==12,1);
%至此EVENT整理完成，[START,target_on,go_signal,move_onset,reach,touch,END,condition,move_onset_two,reach_two,touch_two]

spikedata = spikedata;  %#ok<*ASGSL,*NODEF> % blackrock文件中的电峰时间序列
for i=1:length(event_full)
    spikedata(spikedata(:,1) > event_full(i,1) & spikedata(:,1) < event_full(i,7),2) = i;%标记每个spike所属的trial
end

spikedata(spikedata(:,2)==0,:)=[];%去掉不属于任何trial的spike
% trial_ind = 1:length(event_full); %unique(spikedata(:,2));%有spike的trail
event = event_full;
spike_processed = zeros(size(event,1),sum(spikedata(:,2) == mode(spikedata(:,2))));

for k=1:length(event)
    spike_processed(k,1:sum(spikedata(:,2) == k)) = spikedata(spikedata(:,2) == k,1);
end

spike_processed(spike_processed == 0) = nan; %至此数据整理完毕，主体为spike_processed,event,trial_ind,event_trial_data
end