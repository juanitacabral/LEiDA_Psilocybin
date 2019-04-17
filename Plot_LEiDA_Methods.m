% function Plot_LEiDA_Methods

%
%  This code makes different types of Figures:
%  A - ILLUSTRATE BOLD HILBERT PHASE
%  B - PLOT PHASES IN CORTEX, IN UNIT CIRCLE, PC matrox and V1
%  C - PLOT ALL EIGENVECTORS FROM ONE SUBJECT
%  D - PLOT THE REPERTOIRE OF CENTROIDS FOR A GIVEN K

%% A 

% LOAD BOLD SIGNAL FROM ONE SESSION
Subject=1;
Condition=1; % 3=baseline, 4=psilo
TR = 3;
seed=2; % Choose one AAL seed to plot
load psilotc_total.mat tc_aal
signal = tc_aal{Subject,Condition};
[N_areas, Tmax]=size(signal);
clear tc_aal 

Order=[1:2:N_areas N_areas:-2:2];

% FILTER SETTINGS              
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                    %highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
signal_filt=zeros(size(signal));
Complex_BOLD=complex(zeros(size(signal)));

% Filter signal and get the Hilbert
for seed=1:N_areas
    signal(seed,:)=detrend(signal(seed,:)-mean(signal(seed,:)));
    signal_filt(seed,:) =filtfilt(bfilt,afilt,signal(seed,:));
    Complex_BOLD(seed,:) = hilbert(signal_filt(seed,:));
end

Theta=angle(Complex_BOLD(seed,1:50));
%Amp=abs(Complex_BOLD(seed,1:50));

% Panel A - FIGURE ILLUSTRTING BOLD HILBERT PHASE
figure
plot3(0:TR:(length(Theta)-1)*TR,sin(Theta),cos(Theta),'k','LineWidth',2);
xlabel('Time (seconds)')
ylabel('Imag')
zlabel('Real')
grid on
ylim([-3 3])
zlim([-3 3])
hold on
P1=[0:TR:(length(Theta)-1)*TR;zeros(size(Theta));zeros(size(Theta))]';
P2=[0:TR:(length(Theta)-1)*TR;sin(Theta);cos(Theta)]';
daspect([10 1 1])
arrow3(P1,P2,'Arrowhead',.5)
arrow3([-20 0 0],[length(Theta)*TR+30 0 0])
xlim([-20 length(Theta)*TR+30])
hold on
plot3(0:TR:(length(Theta)-1)*TR,3*ones(size(Theta)),signal_filt(seed,1:50),'b');
plot3(0:TR:(length(Theta)-1)*TR,3*ones(size(Theta)),signal(seed,1:50),'g');
plot3(0:TR:(length(Theta)-1)*TR,3*ones(size(Theta)),cos(Theta),'k:');
plot3(0:TR:(length(Theta)-1)*TR,sin(Theta),-3*ones(size(Theta)),'k:');
P1(:,1)=-20;
P2(:,1)=-20;
arrow3(P1,P2,'Arrowhead',0)

%% Panel B - PLOT CORTEX WITH ARROWS, PHASE PROTRAIT, MATRIX

ThetaB=angle(Complex_BOLD);
clear signal signal_filt Theta P1 P2 Complex_BOLD

for t=1:Tmax
    %Calculate the Instantaneous BOLD Phase Coherence)
    iFC=zeros(N_areas);
    for n=1:N_areas
        for p=1:N_areas
            iFC(n,p)=cos(ThetaB(n,t)-ThetaB(p,t));
        end
    end
    % Get the leading eigenvector
    [vec1(:,t) val1(:,t)]=eigs(iFC,1);
end
vec1(:,sum(vec1)>0)=-vec1(:,sum(vec1)>0);


figure
subplot(2,4,[1 2 5 6])
hold on

% % PLOT CORTEX (THIS PART NEEDS SPM, otherwise comment)
cortex.path='MNI152_T1_2mm_brain_mask.nii';
cortex.pial=mapPial(cortex.path);
cortex.color=[0.9 0.9 0.9];
cortex.transparency=0.3; % To view only opaque cortex =1;
cortex.val=0.2;
redux=1;
sregion=smooth3(cortex.pial);
psregion=patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', cortex.color, 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', cortex.transparency); %transparency

% center origin

ori=[65 45.5 35];
load aal_cog.txt aal_cog
scale=5.5;
MNI_coord=scale*(aal_cog/10);
clear aal_cog

% PLOT PHASES ON CORTEX
t=2; 
V=vec1(:,t);
scale=6;
P1=[MNI_coord(:,2)+ori(1) MNI_coord(:,1)+ori(2) MNI_coord(:,3)+ori(3)];
P2=P1+[scale*cos(ThetaB(:,t)-(pi/2)) scale*sin(ThetaB(:,t)-(pi/2)) zeros(N_areas,1)];
daspect([1 1 1])
blueharrow3=arrow3(P1(V<0,:),P2(V<0,:),'b',.5,.7,.5);
redharrow3=arrow3(P1(V>=0,:),P2(V>=0,:),'r',.5,.7,.5);
axis off;
axis equal
set(gcf,'Renderer', 'OpenGL') % USE UNDER LINUX FOR TRANSPARENCY
view(-90,90)
set(gca,'CameraViewAngle', 6);
set(gca, 'Projection', 'orthographic')
set(gca, 'CameraTarget', [51 68 90])
material dull; lighting phong;
camlight;
rotate3d;

% PLOT PHASES ON COMPLEX PLANE
subplot(2,4,3)
hold on
B1=ThetaB(:,t);
B0=zeros(N_areas,1);
plot([B0(V<0) cos(B1(V<0))]',[B0(V<0) sin(B1(V<0))]','b')
plot([B0(V>=0) cos(B1(V>=0))]',[B0(V>=0) sin(B1(V>=0))]','r')
ylabel('Im')
xlabel('Re')
xlim([-1 1])
ylim([-1 1])
axis square
title('BOLD Phases, \Theta_n(t)')

% PLOT PC MATRIX
subplot(2,4,7)  
iFC=zeros(N_areas);
for n=1:N_areas
    for p=1:N_areas
        iFC(n,p)=cos(ThetaB(n,t)-ThetaB(p,t));
    end
end
colormap(jet)
imagesc(iFC(Order,Order))
ylabel('Brain Areas')
xlabel('Brain Areas')
set(gca,'XTick',10:20:90)
set(gca,'YTick',10:20:90)
xlim([.5 90.5])
ylim([.5 90.5])
axis square
title({'Phase Coherence','dFC(t)'})


subplot(2,4,[4 8])  % EIGENVECTOR
cla
Vo=V(Order);
hold on
barh(Vo.*(Vo<0),'FaceColor',[0  0  1],'EdgeColor','none','Barwidth',.5)
barh(Vo.*(Vo>=0),'FaceColor',[1 0 0],'EdgeColor','none','Barwidth',.5)
ylim([0 91])
xlim([-.15 .15])
set(gca,'YTick',1:90)
load AAL_labels label90
set(gca,'YTickLabel',label90(Order,:),'Fontsize',6)
title({'Leading','Eigenvector','V_1(t)'},'Fontsize',12)
ylabel('Brain Areas','Fontsize',12)

% TO MAKE THE PLOT EVOLVE OVER TIME

for t=3:1:size(ThetaB,2)
    
    subplot(2,4,[1 2 5 6])
    V=vec1(:,t);
    P2=P1+[scale*cos(ThetaB(:,t)-(pi/2)) scale*sin(ThetaB(:,t)-(pi/2)) zeros(N_areas,1)];
    delete(blueharrow3)
    delete(redharrow3)
    blueharrow3=arrow3(P1(V<0,:),P2(V<0,:),'b',.5,1,.5);
    redharrow3=arrow3(P1(V>=0,:),P2(V>=0,:),'r',.5,1,.5);
    
    subplot(2,4,3)
    cla
    hold on
    B1=ThetaB(:,t);
    plot([B0(V<0) cos(B1(V<0))]',[B0(V<0) sin(B1(V<0))]','b')
    plot([B0(V>=0) cos(B1(V>=0))]',[B0(V>=0) sin(B1(V>=0))]','r')
    ylabel('Imag')
    xlabel('Real')
    xlim([-1 1])
    ylim([-1 1])
    axis square
    title('BOLD Phases, \Theta_n(t)')
    
    subplot(2,4,7)
    cla
    for n=1:N_areas
        for p=1:N_areas
            iFC(n,p)=cos(ThetaB(n,t)-ThetaB(p,t));
        end
    end
    colormap(jet)
    imagesc(iFC(Order,Order))
    ylabel('Brain Areas')
    xlabel('Brain Areas')
    set(gca,'XTick',10:20:90)
    set(gca,'YTick',10:20:90)
    xlim([.5 90.5])
    ylim([.5 90.5])
    axis square
    title({'Phase Coherence','dFC(t)'})
    
    subplot(2,4,[4 8])
    cla
    Vo=V(Order);
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0  0  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 0 0],'EdgeColor','none','Barwidth',.5)
    ylim([0 91])
    xlim([-.15 .15])
    set(gca,'YTick',1:90)
    load AAL_labels label90
    set(gca,'YTickLabel',label90(Order,:),'Fontsize',6)
    title({'Leading','Eigenvector','V_1(t)'},'Fontsize',12)
    ylabel('Brain Areas','Fontsize',12)
    
    pause(0.5)
end

%% C 
% CODE TO PLOT THE EIGENVECTORS OF ALL SUBJECTS OVER TIME (TO SELECT SOME TO PLOT)

N=90;
Order=[1:2:N N:-2:2];
%vec=zeros(1,N);
load /Users/joana/Documents/Work/CarolineEric/NewData/LEIDA_results_sept2018.mat X
Tmax=210;
s=1; % Subject
u=1;
figure('Name',['Subject ' num2str(s)])
for t=1:150
    Vo=X((s-1)*Tmax+t,Order);
    subplot(3,50,u)
    u=u+1;
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 N+1])
    xlim([-.15 .15])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis off
end
    
%% D - PLOT CLUSTER CENTROIDS FOR A GIVEN SOLUTION

load LEiDA_psilo_newkresults.mat Kmeans_results
k=4;
figure
for c=1:k
    Vo=Kmeans_results{k}.C(c,Order);
    subplot(1,k,c)
    u=u+1;
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 91])
    xlim([-.15 .15])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis on
end
    

