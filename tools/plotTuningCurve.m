function fig = plotTuningCurve

load('data/exampleFit.mat');
n = 60;

n_sp = squeeze(co_spm(n,:,:));
n_dsp = squeeze(dr_spm(n,:,:));z=[];
err = [squeeze(2*SE_dr_N(n,:,:));squeeze(2*SE_co_N(n,:,:))]';

z = [n_dsp;n_sp]';

%% neural data tuning
T_bin = [43,53];
for p = 1:2
    fig = figure;
    col = [245, 124, 0; 162, 20, 47;183, 70, 255; 0, 114, 189; 119, 172, 48]/255;
    col = col([2,4,3,1,5],:);

    s(1)= subplot(1,3,1);
    for i = 1:5
        errorbar([0:60:360],[z(T_bin(p),(i-1)*6+1:i*6),z(T_bin(p),(i-1)*6+1)],[err(T_bin(p),(i-1)*6+1:i*6),err(T_bin(p),(i-1)*6+1)],'o','markerfacecolor',col(i,:),'markeredgecolor',col(i,:),'LineWidth',2,'color',col(i,:));
        hold on
        [fitresult,~] = createFit([dr_am(1:6),360], [z(T_bin(p),(i-1)*6+1:i*6),z(T_bin(p),(i-1)*6+1)]);
        a = plot( fitresult);set(a,'linewidth',2,'color',col(i,:))
    end
    title('Neuron data');xlim([0,360]);ylim([0,1]);box off;xlabel('First Direction (\circ)');ylabel('Normalized FR')
    axis square;legend off; set(gca,'xtick',0:60:360);
    
    
    
    %% model tuning

    FD = 1:360;
    ccw120SD = FD-210;
    cw120SD = FD-150;
    SD180 = FD-180;
    ccw60SD = FD-240;
    cw60SD = FD-120;
    
    % n = 150;
    Amodel = d_fitmodel_add(n,:);Afr=[];
    Mmodel = d_fitmodel_multi(n,:);Mfr=[];
    for t = 1:91
        Afitmodel = Amodel{t};Mfitmodel = Mmodel{t};
        Afr(t,:) = Afitmodel(repmat(FD,[1,5]), [ccw120SD,cw120SD,SD180,ccw60SD,cw60SD]);
        Mfr(t,:) = Mfitmodel([ccw120SD,cw120SD,SD180,ccw60SD,cw60SD],repmat(FD,[1,5]));
    end
    
    
    % addition
    s(2)= subplot(1,3,2);
    for i = 1:5
        plot(FD,Afr(T_bin(p),(i-1)*length(FD)+1:i*length(FD)),'color',col(i,:),'linewidth',2);hold on
    end
    set(gca,'xtick',0:60:360);
    title('Additive');xlim([0,360]);ylim([0,1]);box off;xlabel('First Direction (\circ)');%ylabel('Normalized FR');
    axis square
    % multi
    s(3)= subplot(1,3,3);
    for i = 1:5
        plot(FD,Mfr(T_bin(p),(i-1)*length(FD)+1:i*length(FD)),'color',col(i,:),'linewidth',2);hold on
    end
    title('Multiplicative');xlim([0,360]);ylim([0,1]);box off;xlabel('First Direction (\circ)');%ylabel('Normalized FR')
    axis square; set(gca,'xtick',0:60:360);
    
    
    set(gcf,'position',[350 350 1300 400]);
end
end

