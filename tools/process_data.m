function [spike_processed,event] = process_data(pathname,filename)

load(strcat(pathname,filename));

%����spta�ã����뼡������

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
reach = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==12 & x(:,2)==0))-1)]*x(find(x(:,3)==12 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);%12��©������
ind2 = cellfun('length',reach)==0;
ind = cellfun('length',move_onset)==0;%�ҳ�reachΪ�յ�trial��ΪCO���ٴ�reach_two�ж�Ӧ��ֵת�Ƶ�reach�У�����ͳһ��CO��reach��DR�ĵ�һ��reach
ind3 = ind2-ind;
correct_rec_beh = correct_rec_beh(ind3==0);
condition = condition(ind3==0);%ȥ��reach©���trial

move_onset_two = cellfun(@(x)x(x(:,3)==13 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
reach_two = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==14 & x(:,2)==0))-1)]*x(find(x(:,3)==14 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
ind2 = cellfun('length',reach_two)==0;
ind = cellfun('length',move_onset_two)==0;%�ҳ�reachΪ�յ�trial��ΪCO���ٴ�reach_two�ж�Ӧ��ֵת�Ƶ�reach�У�����ͳһ��CO��reach��DR�ĵ�һ��reach
ind3 = ind2-ind;
correct_rec_beh = correct_rec_beh(ind3==0);
condition = condition(ind3==0);%ȥ��reach_two©���trial

move_onset = cellfun(@(x)x(x(:,3)==11 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
move_onset_two = cellfun(@(x)x(x(:,3)==13 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
reach = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==12 & x(:,2)==0))-1)]*x(find(x(:,3)==12 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
reach_two = cellfun(@(x)[1,zeros(1,length(find(x(:,3)==14 & x(:,2)==0))-1)]*x(find(x(:,3)==14 & x(:,2)==0)),correct_rec_beh,'UniformOutput',0);
touch = cellfun(@(x)x(x(:,3)==3 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
touch_two = cellfun(@(x)x(x(:,3)==7 & x(:,2)==0),correct_rec_beh,'UniformOutput',0);
ind = cellfun('length',reach)==0;%�ҳ�reachΪ�յ�trial��ΪCO���ٴ�reach_two�ж�Ӧ��ֵת�Ƶ�reach�У�����ͳһ��CO��reach��DR�ĵ�һ��reach
move_onset(ind) = move_onset_two(ind);
move_onset_two(ind) = [];
reach(ind) = reach_two(ind);
reach_two(ind) = [];
touch(ind) = touch_two(ind);
touch_two(ind) = [];
% a=[1,zeros(1,length(find(correct_rec_beh{1,2}(:,3)==14 & correct_rec_beh{1,2}(:,2)==0))-1)]
% b=correct_rec_beh{1,2}(find(correct_rec_beh{1,2}(:,3)==14 & correct_rec_beh{1,2}(:,2)==0),:)
% a*b
% c=cell2mat(c) ����һ����ȡ0 14 0,0 12 0markerΨһ����˼·
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

for i=1:length(spikeTime)%��һ����spike��ʱ�̣��ڶ����Ƕ�Ӧ��trial������������condition��������ȡbin marker��Ϣ
    for k=1:length(event_trial_data)%��Ŀ����һ����
    if event_trial_data(k,1) < spikeTime(i,1) & event_trial_data(k,7) > spikeTime(i,1)
        if event_trial_data(k,2)+mBin(1) < spikeTime(i,1) & event_trial_data(k,2)+nBin(1) > spikeTime(i,1)
            spikeTime(i,4)=1;%��4�С�1����ʾtarget_on
        elseif event_trial_data(k,3)+mBin(2) < spikeTime(i,1) & event_trial_data(k,3)+nBin(2) > spikeTime(i,1)
            spikeTime(i,5)=2;%��5�С�2����ʾgo_signal
        elseif event_trial_data(k,4)+mBin(3) < spikeTime(i,1) & event_trial_data(k,4)+nBin(3) > spikeTime(i,1)
            spikeTime(i,6)=3;%��6�С�3����ʾmove_onset 
        elseif event_trial_data(k,4)+mBin(4) < spikeTime(i,1) & event_trial_data(k,4)+nBin(4) > spikeTime(i,1)
            spikeTime(i,7)=4;%��7�С�4����ʾperimovement
        elseif event_trial_data(k,5)+mBin(5) < spikeTime(i,1) & event_trial_data(k,5)+nBin(5) > spikeTime(i,1)
            spikeTime(i,8)=5;%��8�С�5����ʾreach
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
spike_processed(spike_processed == 0) = nan; %��������������ϣ�����Ϊspike_processed,event,trial_ind,event_trial_data,DR
end

