
% Calculate overlap with the 7 resting-state Networks from Yeo et al. 2011

% load the masks (AAL and Yeo networks) in MNI space 2mm
load('/Users/joana/Documents/Work/LEiDA general/Networks/yeo_aal_2mm.mat','aal2mm','yeo2mm')

N_areas=90;
Yeo_AAL=zeros(7,N_areas);

% Create 7 vectors representing the 7 RSNs in AAL space 
for n=1:90
    indn=find(aal2mm==n);
    for Net=1:7
        Yeo_AAL(Net,n)=numel(find(yeo2mm(indn)==Net));
    end
end

%Figure of the Yeo Repertoire in AAL space

figure
Order=[1:2:N_areas N_areas:-2:2];
for Net=1:7
    subplot(3,7,[Net Net+7])
    hold on
    barh(squeeze(Yeo_AAL(Net,Order)),'EdgeColor','none','Barwidth',.5,'FaceColor','r')
    ylim([0 91])
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',8)    
    set(gca,'YTickLabel',[])
    subplot(3,7,[Net+7*2])
    plot_nodes_in_cortex_yeo(Yeo_AAL(Net,:))
end

load('LEiDA_psilo_newkresults.mat', 'Kmeans_results')
load /Users/joana/Documents/Work/Mushrooms/LEIDA_clouds_files/LEiDA_psilo_data.mat Time_sessions
k=7;
[~, ind_sort]=sort(hist(Kmeans_results{k}.IDX(Time_sessions(1,:)==1),1:k),'descend');
VLeida=Kmeans_results{k}.C(ind_sort,:);

overlap_yeo_aal=zeros(k,7);
overlap_yeo_aal_p=zeros(k,7);

for FLeida=1:k
    V=VLeida(FLeida,:);
    %V=round(V*1e3)/1e3;  % round very small values close to zero
    V=V/max(abs(V));
   
    for NetYeo=1:7
        [cc p]=corrcoef(V.*(V>0),Yeo_AAL(NetYeo,:));
        overlap_yeo_aal(FLeida,NetYeo)=cc(2);
        overlap_yeo_aal_p(FLeida,NetYeo)=p(2);
    end
end

% Figure of the overlap bars

figure
bar(1:k,overlap_yeo_aal)
overlap_yeo_aal=overlap_yeo_aal(2:k,:);
overlap_yeo_aal_p=overlap_yeo_aal_p(2:k,:);

