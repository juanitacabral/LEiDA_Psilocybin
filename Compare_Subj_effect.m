load LEiDA_psilo_data.mat Time_sessions
load('LEiDA_psilo_newkresults.mat')
k=7;
IDX_Psi=Kmeans_results{k}.IDX;

condition=1;
for c=1:k
    ProbC(c)=mean(IDX_Psi(Time_sessions(1,:)==condition)==c);
end
[~, ind_sort]=sort(ProbC,'descend'); 


load ClusterpsiloStatsNew_Paired.mat
P_psi=squeeze(P(6,:,2,ind_sort));

load ClusterplaceboStats.mat
P_pla=squeeze(P(:,2,ind_sort));

Probs_Psi_Pla=cat(2,P_psi',P_pla');

load Effect_subjs Effects_Psi_Pla

corrPC_Effect=zeros(1,7);
pval=zeros(1,7);

k=7;
for c=1:k
    if find([2 4 5 6 7]==c) % For these networks we have no a priori hypothesis
    [corrPC_Effect(c), pval(c)]=corr(squeeze(Probs_Psi_Pla(c,:))',Effects_Psi_Pla','Tail','both');
    elseif c==1   % Because we know this state increases with psilocybin
        [corrPC_Effect(c), pval(c)]=corr(squeeze(Probs_Psi_Pla(c,:))',Effects_Psi_Pla','Tail','right');
    elseif c==3   % Because we know this state decreases with Psilocybin
        [corrPC_Effect(c), pval(c)]=corr(squeeze(Probs_Psi_Pla(c,:))',Effects_Psi_Pla','Tail','left');
    end
end

figure
for c=1:7
    subplot(1,7,c)
    hold on
    plot(Probs_Psi_Pla(c,1:9),Effects_Psi_Pla(1:9),'r*')
    plot(Probs_Psi_Pla(c,10:18),Effects_Psi_Pla(10:18),'b*')
    title({['cc=' num2str(round(corrPC_Effect(c)*1000)/1000)],['p=' num2str(round(pval(c)*10000)/10000)]})
    for s=1:9
        plot([Probs_Psi_Pla(c,s) Probs_Psi_Pla(c,s+9)],[Effects_Psi_Pla(s) Effects_Psi_Pla(s+9)],'k')
    end
    xlabel('PC state Prob.')
    if c==1
        ylabel('Subjective Rating')
    end
end


% corrPC_Effect =
% 
%     0.1699    0.3555   -0.5558   -0.0817   -0.1291   -0.3749    0.2135
%     
% pval =
% 
%     0.5003    0.1476    0.0166    0.7474    0.6096    0.1252    0.3950
