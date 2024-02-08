% P test
load('D:\OneDrive\Documents\PaperMeta\code\analysis_by_WTW\NN_Rvision\perm_dir\Perm_result_fitnlm_1st.mat')
perm_rs(:,:,:,1) = d_rs_full;
perm_coff(:,:,:,:,1) = d_coff_full;
load('D:\OneDrive\Documents\PaperMeta\code\analysis_by_WTW\NN_Rvision\perm_dir\Perm_result_fitnlm_mul.mat')
perm_rs(:,:,:,2) = d_rs_full;
perm_coff(:,:,:,:,2) = d_coff_full;
load('D:\OneDrive\Documents\PaperMeta\code\analysis_by_WTW\NN_Rvision\perm_dir\Perm_result_fitnlm_2nd.mat')
perm_rs(:,:,:,3) = d_rs_full;
perm_coff(:,:,:,:,3) = d_coff_full;

clearvars -except perm_rs perm_coff t

rs_m = squeeze(nanmean(perm_rs,3));
coff_m = squeeze(mean(abs(perm_coff),4));


% M1
for p = 1:3
    rs = rs_m(41:end,:,p);
    coff = coff_m(41:end,:,:,p);
    
    %Coeff
    col=[21,112,177;29,153,29;255,119,0;118 113 113;141,194,211]/255;
    a = squeeze(mean(abs(coff(:,:,[1 3 4 5])),1));
    coff = abs(coff(:,:,[1 3 4 5]));E=[];
    for i=1:4
        E(:,i) = squeeze(std(coff(:,:,i),1)/sqrt(size(coff(:,:,i),1)));
    end
    
    s(6)=subplot(4,3,6);
    fill([t fliplr(t)],[a(:,p)'+2*E(:,p)' fliplr(a(:,p)'-2*E(:,p)')],col(p,:),'LineStyle','none');alpha(.5);
    hold on
    
end
s(6).Position(2)=0.58;

% pmD
for p = 1:3
    rs = rs_m(1:40,:,p);
    coff = coff_m(1:40,:,:,p);
    
    %Coeff
    col=[21,112,177;29,153,29;255,119,0;118 113 113;141,194,211]/255;
    a = squeeze(mean(abs(coff(:,:,[1 3 4 5])),1));
    coff = abs(coff(:,:,[1 3 4 5]));E=[];
    for i=1:4
        E(:,i) = squeeze(std(coff(:,:,i),1)/sqrt(size(coff(:,:,i),1)));
    end
    
    s(12)=subplot(4,3,12);
    fill([t fliplr(t)],[a(:,p)'+2*E(:,p)' fliplr(a(:,p)'-2*E(:,p)')],col(p,:),'LineStyle','none');alpha(.5);
    hold on
    
end
