
load ClusterpsiloStatsNew_Paired.mat pvalue_LvsP P
rangeK_run=2:20;
rangeK=5:10;

P_pval=pvalue_LvsP(find(rangeK_run==rangeK(1)):find(rangeK_run==rangeK(end)),1:rangeK(end));

sigC=zeros(1,length(rangeK));
Min_p_value=zeros(1,length(rangeK));

for k=1:length(rangeK)   
    [Min_p_value(k), sigC(k)]=min(P_pval(k,P_pval(k,:)>0));    
end

figure
semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)
        if P_pval(k,c)<0.05 && P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),P_pval(k,c),'*r');
        end
        if P_pval(k,c)>0.05
            semilogy(rangeK(k),P_pval(k,c),'*k');
        end
        if P_pval(k,c)<(0.05/rangeK(k)) && P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*g');
        end
        if P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*b');
        end
    end
end

ylabel('Prob before vs after psilocybin (p-value)')
xlabel('Number of states k')
set(gca,'XTick', rangeK)
ylim([1e-4 1])
xlim([rangeK(1)-1 rangeK(end)+1])
box off

load LEiDA_psilo_newkresults.mat Kmeans_results

figure % FIGURE SHOWING MOST DIFFERENT NETWORK IN EACH K
load AAL_labels.mat label90
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

V_all=zeros(N_areas,length(rangeK));

for k=1:length(rangeK)
    c=sigC(k);
    V=Kmeans_results{rangeK(k)}.C(c,:);
    V=V(Order)/max(abs(V));
    V_all(:,k)=V;
    
    subplot(1,length(rangeK),k)
    hold on
    barh(V.*(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(V.*(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 91])
    xlim([-1 1])
    
    %set(gca,'Ydir','reverse')
    
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',8)    
%     if k==1
%         set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
%     else
        set(gca,'YTickLabel',[])
%     end

    if Min_p_value(k)<0.05 && Min_p_value(k)>(0.05/rangeK(k))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','r')
    elseif Min_p_value(k)>0.05
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','k')
    elseif Min_p_value(k)<(0.05/k) 
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','g')
    end

end

cc_all=corrcoef(V_all)
% 
% 
% 
% 
% 
% 
% figure
% semilogy(rangeK,0.05*ones(1,length(rangeK)),'r--')
% hold on
% semilogy(rangeK,0.05./rangeK,'g--')
% 
% for k=1:length(rangeK)   
%     semilogy(rangeK(k),P_Prob(k,P_Prob(k,:)>0),'*k');    
% end
% semilogy(rangeK,Min_p_value,'*b')
% %semilogy(rangeK,Min_p_value.*rangeK,'g*-')
% semilogy(5:17,0.05./sum(5:17)*ones(1,length(5:17)),'b--')
% 
% semilogy(2:20,0.05./sum(2:20)*ones(1,length(2:20)),'c--')
% 
% ylabel('Prob Patients vs Controls (p-value)')
% xlabel('Number of clusters K')
% box off
% 
% % LIFETIMES
% 
% sigC=zeros(length(rangeK));
% 
% for k=1:length(rangeK)   
%     [Min_p_value(k), sigC(k)]=min(P_LT(k,P_LT(k,:)>0));    
% end
% 
% figure
% semilogy(rangeK,0.05*ones(1,length(rangeK)),'r--')
% hold on
% semilogy(rangeK,0.05./rangeK,'g--')
% 
% for k=1:length(rangeK)   
%     semilogy(rangeK(k),P_LT(k,P_LT(k,:)>0),'*k');    
% end
% semilogy(rangeK,Min_p_value,'*b')
% %semilogy(rangeK,Min_p_value.*rangeK,'g*-')
% 
% ylabel('LT Patients vs Controls (p-value)')
% xlabel('Number of clusters K')
% box off
% 
% figure
% load /Users/joana/Documents/Work/CarolineEric/CodesData/AAL_labels.mat label90
% load /Users/joana/Documents/Work/CarolineEric/NewData/LEiDA_results_sept2018.mat Kmeans_results
% N_areas=90;
% Order=[1:2:N_areas N_areas:-2:2];
% label90=label90(Order,:);
% 
% V_all=zeros(N_areas,length(rangeK));
% 
% for k=1:length(rangeK)
%     c=sigC(k);
%     V=Kmeans_results{rangeK(k)}.C(c,:);
%     V=V(Order);
%     V_all(:,k)=V;
%     
%     subplot(1,length(rangeK),k)
%     hold on
%     barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
%     barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
%     ylim([0 91])
%     xlim([-.15 .15])
%     
%     grid on
%     set(gca,'YTick',1:N_areas,'Fontsize',8)    
%     if k==1
%         set(gca,'YTickLabel',label90(end:-1:1,:))
%     else
%         set(gca,'YTickLabel',[])
%     end
%     title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10)
% end
% 
% V_mean=mean(V_all(:,Min_p_value<0.05),2)';
% 
% figure
% subplot(3,2,1)
% plot_nodes_in_cortex(V_mean)
% rotate3d
% view(0,0)
% 
% subplot(3,2,3)
% plot_nodes_in_cortex(V_mean)
% view(-90,0)
% rotate3d
% 
% subplot(3,2,[2 4 6])
% hold on
% barh(find(V_mean<0),V_mean(V_mean<0),'FaceColor',[.5 .5 .5],'EdgeColor','none')
% barh(find(V_mean>=0),V_mean(V_mean>=0),'FaceColor','red','EdgeColor','none','Barwidth',.4)
% set(gca,'YTick',1:90,'Fontsize',7)
% set(gca,'YTickLabel',label90(end:-1:1,:))
% ylim([0 91])
% xlim([-.15 .15])
% title('relative BOLD phase','Fontsize',12)
% grid on
% 
% subplot(3,2,5)
% colormap(jet)
% FC_V=V_mean'*V_mean;  
% li=max(abs(FC_V(:)));
% imagesc(FC_V,[-li li])
% axis square
% title('FC pattern') 
% ylabel('Brain area #')
% xlabel('Brain area #')  