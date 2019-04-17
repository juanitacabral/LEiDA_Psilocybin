


load LEiDA_psilo_data.mat Time_sessions
load LEiDA_psilo_newkresults.mat

% For each K - Plot the PC state as network in cortex
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];

u=0;
figure
for K=5:10
    k=K;
    % Get the K patterns
    
    [~, ind_sort]=sort(hist(Kmeans_results{k}.IDX(Time_sessions(1,:)==1),1:k),'descend');
    V=Kmeans_results{k}.C(ind_sort,:);
    
    
    for c=1:K
        
        subplot(6,10,c+u*10)
        % This needs function plot_nodes_in_cortex.m and aal_cog.m
        plot_nodes_in_cortex(V(c,:))
    end
    u=u+1;
end
