function LEiDA_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function processes the data for LEiDA
%
% - Reads the BOLD data from the folders
% - Computes the BOLD phases
% - Calculates the Order Parameter, Mestastability and Synchrony
% - Calculates the instantaneous BOLD synchronization matrix
% - Calculates the instantaneous Leading Eigenvector
% - Calculates de FCD matrices
%
%  The data is from 9 healthy subjects
%    placebo before injection (pcb01)
%    placebo after injectio (pcb02)
%    psilocybin before injection (psi01)
%    psilocybin after injection (psi02)
%
% Saves the Leading_Eigenvectors and FCD matrices into LEiDA_data.mat
%
% Joana Cabral May 2016 
% joana.cabral@psych.ox.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load psilotc_total.mat tc_aal

[n_Subjects, n_Task]=size(tc_aal);
[N_areas, Tmax]=size(tc_aal{1,1});


% Create empty variable to save patterns 
Leading_Eig=cell(9,4);
FCD_eig=cell(9,4);
FCD_iFC=cell(9,4);
iFC_values=zeros(Tmax,(N_areas*(N_areas-1)/2));

% FILTER SETTINGS
TR=2;                 
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                   %highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

for s=1:n_Subjects
    
    for task=1:n_Task
        
        signal = tc_aal{s,task};
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
            ts=demean(detrend(signal(seed,:)));
            signal_filt =filtfilt(bfilt,afilt,ts);
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
        
        %Calculate the Kuramoto Order Parameter
        OP=abs(sum(exp(1i*Phase_BOLD))/N_areas);
        Metasta(s,task)=std(OP);
        Synchro(s,task)=mean(OP);
        
        for t=1:Tmax
            
            %Calculate the Instantaneous FC (BOLD Phase Coherence)
            iFC=zeros(N_areas);
            
            for n=1:N_areas
                for p=1:N_areas
                    iFC(n,p)=cos(adif(Phase_BOLD(n,t),Phase_BOLD(p,t)));
                end
            end
            
            
            
            
            
            
            
            % Get the leading eigenvector
            [Leading_Eig{s,task}(t,:),~]=eigs(iFC,1);
            iFC_values(t,:)=iFC(triu(ones(N_areas),1)>0);
        end
        
        % Calculate the FCD
        for t=1:Tmax
            eig1=squeeze(Leading_Eig{s,task}(t,:));
            iFC1=iFC_values(t,:);
            for t2=1:Tmax
                eig2=squeeze(Leading_Eig{s,task}(t2,:));
                iFC2=iFC_values(t2,:);
                % Cosine similarity between vectors at t1 and t2
                FCD_eig{s,task}(t,t2)=dot(eig1,eig2)/norm(eig1)/norm(eig2);
                FCD_iFC{s,task}(t,t2)=dot(iFC1,iFC2)/norm(iFC1)/norm(iFC2);
            end
        end
    end
end

save('LEiDA_data','FCD_eig','FCD_iFC','Metasta','Synchro','Leading_Eig')




