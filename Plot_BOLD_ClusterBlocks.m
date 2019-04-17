%% script to plot the blocks of LEiDA states over BOLD

load LEiDA_Task_kmeans.mat Kmeans_results
load LEiDA_Task_data.mat Time_subjects

% Nr of states
K=8;
k=8;

Directory='/Users/joana/Documents/Work/Eloise/';
names=dir([Directory 'BRST_AAL90_tc/*.txt']);

N_areas = 90;

s=1;
BOLD=load([Directory 'BRST_AAL90_tc/' names(s).name]);

T=Time_all==s;
Ctime=Kmeans_results{k}.IDX(T);

Tmax=280-2;
TR=2.2;

Directory='/Users/joana/Documents/Work/Eloise/timings/';
timing_folders=dir([Directory '*_timing']);
n_Subjects=length(timing_folders);

Baseline_start=zeros(1,n_Subjects);

N_TR_after_task=2;


% BAND-PASS FILTER SETTINGS            
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter (Hz)
fhi = 0.1;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

        for seed=1:N_areas
            ts=BOLD(seed,:)-mean(BOLD(seed,:));
            BOLD(seed,:) =filtfilt(bfilt,afilt,ts);
        end

timing_all=zeros(1,n_Subjects*Tmax);
for s=1:n_Subjects
    timings=load([Directory timing_folders(s).name '/ev1.txt']);
    % Baseline starts at last event + 3 TR (set as -1)
    Baseline_start(s)=round((timings(end,1)+timings(end,2))/TR)-1+3;
    timing_all((s-1)*Tmax + Baseline_start(s) : (s-1)*Tmax + Tmax)=-1; 
    
    for event=1:6
        
    % Timings of the events
    % ev1 = 1, ev2=2,..., between events marked as 0)
    %    t_events=round((timings(:,1))/TR);
    timings=load([Directory timing_folders(s).name '/ev' num2str(event) '.txt']);
    timing_all((s-1)*Tmax + round(timings((timings(:,3)==1),1)/TR)+N_TR_after_task-1)=event;
    
    end
end






    
%     % FILTER SETTINGS
%     TR=2;
%     fnq=1/(2*TR);                 % Nyquist frequency
%     flp = .02;                    % lowpass frequency of filter
%     fhi = 0.1;                    %highpass
%     Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
%     k=2;                          % 2nd order butterworth filter
%     [bfilt,afilt]=butter(k,Wn);
     ROISignals=detrend(ROISignals-mean(ROISignals));
     signal = ROISignals';
     taskln=Tmax; % length of current condition
%     signal_filt= [];
%     
%     % Filter the BOLD
%     for seed=1:N_areas
%         ts=signal(seed,:)-mean(signal(seed,:));
%         signal_filt(seed,:) =filtfilt(bfilt,afilt,ts);
%     end
    
    % Create Figure BOLD with shaded State Activations on white
    % background
    
    figure('color','w')
    
    % first plot the states as a stairplot
    subplot(2,1,2)
    stairs(1+0.8*(Ctime==1),'b','LineWidth',2)
    hold on
    stairs(2+0.8*(Ctime==2),'r','LineWidth',2)
    stairs(3+0.8*(Ctime==3),'g','LineWidth',2)
    set(gca,'YTick',1:3)
    set(gca,'YTickLabel',['State 1';'State 2';'State 3'],'Fontsize',12)
    ylim([0.5 4])
    xlim([0 taskln])
    
        % then plot the states
        subplot(2,1,1)
        hold on
        ymax=150;
        xlim([0 taskln])
        
        % Blue Blocks
        NetBlue=Ctime==1;
        if Ctime(1)==1 % if this event is on at beginning
            NetON=1;
            NetOF_final=taskln;
        else
            NetON=[];
            NetOF_final=[];
        end
        NetON=[NetON; find(diff(NetBlue)==1)'];
        NetOF=[find(diff(NetBlue)==-1)'; NetOF_final];
        if size(NetOF,1)<size(NetON,1) % balance for final event
            NetOF=[NetOF; taskln];
        end
        for t=1:length(NetON)
            ON=NetON(t);
            if NetOF(t)>ON
                OFF=NetOF(t);
            else
                if t==length(NetON)
                    OFF=taskln;
                else
                    OFF=NetOF(t+1);
                end
            end
            x = [ON OFF OFF ON];
            y = [-ymax -ymax ymax ymax];
            p=patch(x,y,'r');
            set(p,'LineStyle','none','FaceColor',[0 0 1],'FaceAlpha',0.2);
        end
        
        % Red Blocks
        NetRed=Ctime==2;
        if Ctime(1)==2 % if this event is on at beginning
            NetON=1;
            NetOF_final=taskln;
        else
            NetON=[];
            NetOF_final=[];
        end
        NetON=[NetON; find(diff(NetRed)==1)'];
        NetOF=[find(diff(NetRed)==-1)'; NetOF_final];
        if size(NetOF,1)<size(NetON,1) % balance for final event
            NetOF=[NetOF; taskln];
        end
        for t=1:length(NetON)
            ON=NetON(t);
            if NetOF(t)>ON
                OFF=NetOF(t);
            else
                if t==length(NetON)
                    OFF=taskln;
                else
                    OFF=NetOF(t+1);
                end
            end
            x = [ON OFF OFF ON];
            y = [-ymax -ymax ymax ymax];
            p=patch(x,y,'r');
            set(p,'LineStyle','none','FaceColor',[1 0 0],'FaceAlpha',0.2);
        end
        
        % Green Rectangles
        NetGreen=Ctime==3;
        if Ctime(1)==3 % if this event is on at beginning
            NetON=1;
            NetOF_final=taskln;
        else
            NetON=[];
            NetOF_final=[];
        end
        NetON=[NetON; find(diff(NetGreen)==1)'];
        NetOF=[find(diff(NetGreen)==-1)'; NetOF_final];
        if size(NetOF,1)<size(NetON,1) % balance for final event
            NetOF=[NetOF; taskln];
        end
        for t=1:length(NetON)
            ON=NetON(t);
            if NetOF(t)>ON
                OFF=NetOF(t);
            else
                if t==length(NetON)
                    OFF=taskln;
                else
                    OFF=NetOF(t+1);
                end
            end
            x = [ON OFF OFF ON];
            y = [-ymax -ymax ymax ymax];
            p=patch(x,y,'r');
            set(p,'LineStyle','none','FaceColor',[0 1 0],'FaceAlpha',0.2);
        end
        
        % Plot the BOLD signals in the end
        plot(BOLD(2:end-1,:)')   
        ylabel('BOLD','Fontsize',12)
        xlabel('Time (TR)','Fontsize',12)
        ylim([-max(abs(BOLD(:))) max(abs(BOLD(:)))])
        xlim([0 280])
        set(gca,'Fontsize',12)
        
end

