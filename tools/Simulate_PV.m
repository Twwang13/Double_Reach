function fig = Simulate_PV
%%% surogate data of multi and additive model and test PV performance

%% PLOS comp. Bio. 2016 Hans method
% tau(n) is neuron response latency (normally distributed)
% sigma duration parameter of movement activity (fixed)
% b(n,c) Gain factor (Directional Tuning coeff)
% theta_c direction of reach target
% pd(n) preferred direction
% phy preparatory magnitufe constant (default 0.2)
% mu constant baseline given by mu=sigma*sqrt(-2*ln(phy))
% noise random noise from normal distribution (default SD=0.01)
rng(2021)
col=[1,0,0;0.4 0.7 0.3;0 0.5 0.8;1 0.7 0;0,0.8,0.8;0.8,0,0.8];
sigma=100;%300ms movement activity
mo=200;
tau=normrnd(mo,50,200,1);
pd=360*rand(200,1);
theta=0:60:300;
theta2=theta-210;
phy=0.2;
mu=sigma*sqrt(-2*log(phy));

%%% SR
f_co=[];
for n=1:200
    for c=1:6
        b(n,c)=(1+cosd(theta(c)-pd(n)))/2;
        for t=1:600
            noise = normrnd(5,0.01);
            if t>=tau(n)
                f_co(c,t,n)=b(n,c)*exp(-((t-tau(n)-mu)^2)/(2*(sigma^2)))+noise;
            else
                f_co(c,t,n)=phy*b(n,c)+noise;
            end
        end
    end
end


%% adapt hans model to multi and add model
% tau(n) is neuron response latency (normally distributed)
% sigma duration parameter of movement activity (fixed)
% b(n,c) Gain factor (Directional Tuning coeff)
% theta_c direction of reach target
% pd(n) preferred direction
% phy preparatory magnitufe constant (default 0.2)
% mu constant baseline given by mu=sigma*sqrt(-2*ln(phy))
% noise random noise from normal distribution (default SD=0.01)

% ADDITIVE MODEL
f_add=[];
for n=1:200
    for c=1:6
        b(n,c)=(1+cosd(theta(c)-pd(n))+cosd(theta2(c)-pd(n)))/3;
        for t=1:600
            noise = normrnd(5,0.01);
            if t>=tau(n)
                f_add(c,t,n)=b(n,c)*exp(-((t-tau(n)-mu)^2)/(2*(sigma^2)))+noise;
            else
                f_add(c,t,n)=phy*b(n,c)+noise;
            end
        end
    end
end


%%% multi model
f_mul=[];
for n=1:200
    for c=1:6
        b(n,c)=(1+cosd(theta(c)-pd(n))+cosd(theta(c)-pd(n))*cosd(theta2(c)-pd(n)))/3;
        for t=1:600
            noise = normrnd(5,0.01);
            if t>=tau(n)
                f_mul(c,t,n)=b(n,c)*exp(-((t-tau(n)-mu)^2)/(2*(sigma^2)))+noise;
            else
                f_mul(c,t,n)=phy*b(n,c)+noise;
            end
        end
    end
end
%% surrogate data present plot
model = {'SR','Additive','Multi'};
Data = {f_co,f_add,f_mul}; 
N=3; % choose the example neuron
for m=1:3
    D = Data{m};
    fig = figure;
    set(gcf,'position',[350 250 900 500]);
    subplot(2,3,1)
    for c=1:6
        plot(D(c,:,N)','color',col(c,:),'LineWidth',1);hold on
    end
    box off;set(gca,'xtick',0:600:600);set(gca,'ytick',[]);ylabel('Conditional FR');
    ylim([4.5 6])
     
    % tunning curve
    subplot(2,3,2)
    theta = 1:360;theta2=theta-210;
    tunC = (1+cosd(theta-pd(N)))/2;
    plot(tunC,'--k','LineWidth',2);hold on
    
    switch m
        case 1
            plot(tunC,'-k','LineWidth',2);hold on
        case 2
            tunA = (1+cosd(theta-pd(N))+cosd(theta2-pd(N)))/3;
            plot(tunA,'-k','LineWidth',2);
        case 3
            tunM = (1+cosd(theta-pd(N))+cosd(theta-pd(N)).*cosd(theta2-pd(N)))/3;
            plot(tunM,'-k','LineWidth',2);
    end
   xlim([0,360]);box off
   legend off; set(gca,'xtick',0:120:360);set(gca,'ytick',[])
   ax=gca;set(gca,'color','none');ylim([-0.5 1.5])
   ax.YAxis.Visible = 'off';   % 设置y轴不可见
    
    subplot(2,3,3)
    aa=reshape(permute(D,[3,2,1]),200,600*6);
    [~,score,~] = pca(aa');
    for c=1:6
        plot(score((c-1)*600+1:c*600,1),score((c-1)*600+1:c*600,2),'color',col(c,:),'linewidth',1.5);hold on
    end
    axis off

    %% PV plot
    
    [~,~,PV]=PV_sum(D,pd,mo);
    
    sum_PV=zeros(1,12);
    for i=1:6
        sum_PV=sum_PV+PV(i,:).*exp(1i*(1-i)*pi/3);
    end
    sum_PV_dir_single=angle(sum_PV)/pi*180;
    sum_PV_r_single=abs(sum_PV);
    max_sum_PV_r=max(sum_PV_r_single);
    
    subplot(2,1,2);
    axis([0 600 -108 100]);
    nn=get(gca);
    
    plot_PV_dir=sum_PV_r_single/max_sum_PV_r*100; %长度
    for i=1:size(sum_PV_r_single,2)
        hold on;
        annotation('arrow',[(i-1)*50 (i-1)*50+plot_PV_dir(i)*cosd(sum_PV_dir_single(i)+90)]*nn.Position(3)/600+nn.Position(1),([0 plot_PV_dir(i)*sind(sum_PV_dir_single(i)+90)]+108)*nn.Position(4)/208+nn.Position(2),'LineWidth',1,'HeadLength',10,'HeadLength',10,'Color','k');
    end
    
    plot([0 600],[0 0],'k:');
    axis off;
    plot([0 600],[-100 -100],'k','LineWidth',2);
    plot([0 0],[-98 -108],'k','LineWidth',4);
    plot([600 600],[-98 -108],'k','LineWidth',4);
    text(-10,-130,'0','FontSize',12);
    text(600,-130,'600','FontSize',12);
    
end
end