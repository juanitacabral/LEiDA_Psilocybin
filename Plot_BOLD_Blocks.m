function Plot_BOLD_Blocks

Subject=8;
Conditions=[1 2; 3 4]; %1=pre-pcb, 2=post-pcb 3=baseline, 4=psilo
TR = 3;

%load data
load psilotc_total.mat tc_aal
[N_areas, Tmax]=size(tc_aal{1,1});

% FILTER SETTINGS              
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                    %highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

total_time=Tmax+60/TR+Tmax;
Sessions_BOLD_filt=zeros(2,N_areas,total_time);

load /Users/joana/Documents/Work/Mushrooms/LEIDA_clouds_files/LEiDA_psilo_data.mat Time_sessions
Time_Sub_Pre=find(((Time_sessions(1,:)==1) + (Time_sessions(2,:)==Subject))==2);
Time_Sub_Post=find(((Time_sessions(1,:)==2) + (Time_sessions(2,:)==Subject))==2);

k=7;

load LEiDA_psilo_newkresults.mat Kmeans_results
IDX_Psi_pre=Kmeans_results{k}.IDX(Time_Sub_Pre);
IDX_Psi_post=Kmeans_results{k}.IDX(Time_Sub_Post);

[~, ind_sort]=sort(hist(Kmeans_results{k}.IDX(Time_sessions(1,:)==1),1:k),'descend');
[~,idx_sort]=sort(ind_sort,'ascend');

IDX_Psi_pre=idx_sort(IDX_Psi_pre);
IDX_Psi_post=idx_sort(IDX_Psi_post);
        
clear Kmeans_results

load LEiDA_placebo_kmeans_correct.mat Kmeans_placebo
IDX_Pla_pre=Kmeans_placebo.IDX(Time_Sub_Pre);
IDX_Pla_post=Kmeans_placebo.IDX(Time_Sub_Post);


IDX_Pla_pre=idx_sort(IDX_Pla_pre);
IDX_Pla_post=idx_sort(IDX_Pla_post);
clear Kmeans_results

cmap=[0 0 1; .7 .7 .7 ; 1 0 0 ; 1 0.5 0; 0 1 1 ; 1 0 1; 1 1 0];

for Sessions=1:2 % 1 placebo, 2 psilocybin
    for Cond=1:2 % 1 before, 2 after
        
        signal = tc_aal{Subject,Conditions(Sessions,Cond)};
        signal_filt=zeros(size(signal));

        for seed=1:N_areas
            signal(seed,:)=detrend(signal(seed,:)-mean(signal(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,signal(seed,:));
        end
    
        if Cond==1
            T=1:Tmax;
        else
            T=Tmax+60/TR+1:total_time;
        end
        
        Sessions_BOLD_filt(Sessions,:,T)=signal_filt;
    end
end
   

% FIGURE
figure

ymax=max(abs(Sessions_BOLD_filt(:)));

subplot(5,1,1) % Psilocybin Session
hold on
for t=2:Tmax-1
    x = [t t+1 t+1 t];
    x=(x-1).*TR;
    y = [-ymax -ymax ymax ymax];
    p=patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',cmap(IDX_Psi_pre(t),:),'FaceAlpha',0.5);
end

for t=2:Tmax-1
    x = [t t+1 t+1 t];
    x=(x-1+Tmax+20).*TR;
    y = [-ymax -ymax ymax ymax];
    p=patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',cmap(IDX_Psi_post(t),:),'FaceAlpha',0.5);
end


% Plot the BOLD signals before injection
plot(1:TR:Tmax*TR,squeeze(Sessions_BOLD_filt(2,:,1:Tmax)),'color',[.5 .5 .5])
hold on
% Plot the BOLD signals before injection
plot(Tmax*TR+60+1:TR:(total_time)*TR,squeeze(Sessions_BOLD_filt(2,:,Tmax+60/TR+1:total_time)),'color',[.5 .5 .5])
ylabel({'BOLD','[0.02-0.1Hz]'},'Fontsize',10)
xlabel('Time (seconds)','Fontsize',10)
xlim([0 total_time*TR])
ylim([-ymax ymax])

title(['Subject ' num2str(Subject)])

subplot(5,1,2) % Placebo Session
hold on
for t=2:Tmax-1
    x = [t t+1 t+1 t];
    x=(x-1).*TR;
    y = [-ymax -ymax ymax ymax];
    p=patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',cmap(IDX_Pla_pre(t),:),'FaceAlpha',0.5);
end

for t=2:Tmax-1
    x = [t t+1 t+1 t];
    x=(x-1+Tmax+20).*TR;
    y = [-ymax -ymax ymax ymax];
    p=patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',cmap(IDX_Pla_post(t),:),'FaceAlpha',0.5);
end

% Plot the BOLD signals before injection
plot(1:TR:Tmax*TR,squeeze(Sessions_BOLD_filt(1,:,1:Tmax)),'color',[.5 .5 .5])
hold on
% Plot the BOLD signals before injection
plot(Tmax*TR+60+1:TR:(total_time)*TR,squeeze(Sessions_BOLD_filt(1,:,Tmax+60/TR+1:total_time)),'color',[.5 .5 .5])
ylabel({'BOLD','[0.02-0.1Hz]'},'Fontsize',10)
xlabel('Time (seconds)','Fontsize',10)
xlim([0 total_time*TR])
ylim([-ymax ymax])
  

load ClusterpsiloStatsNew_Paired.mat P
P_Psi=squeeze(P(6,:,:,1:7));

load ClusterplaceboStats P
P_Pla=P;


subplot(5,2,7) % Probabilities before Psilocybin
    P_psi1=squeeze(P_Psi(:,1,ind_sort));
    hold on    
    for c=1:7
        bar((c-0.4)+(1:9)./12,P_psi1(:,c),'FaceColor',cmap(c,:));
    end
    title('Rest before Psilocybin injection')
    ylim([0 0.7])
    ylabel('Probability')
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',{'State 1','State 2','State 3','State 4','State 5','State 6','State 7'},'Fontsize',10)
    
 subplot(5,2,8) % Probabilities under Psilocybin
    P_psi2=squeeze(P_Psi(:,2,ind_sort));
    hold on    
    for c=1:7
        bar((c-0.4)+(1:9)./12,P_psi2(:,c),'FaceColor',cmap(c,:));
    end
    title('Under Psilocybin effects')
    ylim([0 0.7])
    ylabel('Probability')
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',{'State 1','State 2','State 3','State 4','State 5','State 6','State 7'},'Fontsize',10)
 
  
 subplot(5,2,9) % Probabilities before Placebo
    P_pla1=squeeze(P_Pla(:,1,ind_sort));
    hold on    
    for c=1:7
        bar((c-0.4)+(1:9)./12,P_pla1(:,c),'FaceColor',cmap(c,:));
    end
    title('Rest before Placebo injection')
    ylim([0 0.7])
    ylabel('Probability')
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',{'State 1','State 2','State 3','State 4','State 5','State 6','State 7'},'Fontsize',10)
 
    
 subplot(5,2,10) % Probabilities under Psilocybin
    P_pla2=squeeze(P_Pla(:,2,ind_sort));
    hold on    
    for c=1:7
        bar((c-0.4)+(1:9)./12,P_pla2(:,c),'FaceColor',cmap(c,:));
    end
    title('Under Placebo effects')
    ylim([0 0.7])
    ylabel('Probability')
    set(gca,'XTick',1:7)
    set(gca,'XTickLabel',{'State 1','State 2','State 3','State 4','State 5','State 6','State 7'},'Fontsize',10)
 