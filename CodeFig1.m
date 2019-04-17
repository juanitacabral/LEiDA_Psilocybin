
K=7;
k=K;
% Get the K patterns
load LEiDA_psilo_data.mat Time_sessions
[~, ind_sort]=sort(hist(Kmeans_results{k}.IDX(Time_sessions(1,:)==1),1:k),'descend');
V=Kmeans_results{k}.C(ind_sort,:); 
N_areas=90;

figure
colormap(jet) 
% Pannel A - Plot the PC state as a bar plot
% Pannel B - Plot the PC state as network in cortex
% Pannel C - Plot the PC state as matrix
Order=[1:2:N_areas N_areas:-2:2];

for c=1:K
    Vc=V(c,Order);
    Vc=Vc/max(abs(Vc));
    subplot(5,K,[c K+c 2*K+c])
    hold on
    barh(Vc.*(Vc<=0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vc.*(Vc>0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 91])
    xlim([-1 1])
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',8)    
    set(gca,'YTickLabel',[])
    title({['PC State ' num2str(c)],['V_C_' num2str(c)]})

    subplot(5,K,3*K+c)       
    % This needs function plot_nodes_in_cortex.m and aal_cog.m
    plot_nodes_in_cortex(V(c,:))
        
    subplot(5,K,4*K+c)  

    FC_V=Vc'*Vc;  
    li=max(abs(FC_V(:)));
    imagesc(FC_V,[-li li])   
    axis square
    title(['V_C_' num2str(c) '.V_C_' num2str(c) '^T ']) 
    ylabel('Brain area #')
    xlabel('Brain area #')   
end
