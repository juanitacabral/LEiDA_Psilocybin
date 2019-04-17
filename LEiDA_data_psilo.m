function LEiDA_data_VarV1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function preprocesses Psilocybin and Placebo datasets for LEiDA
%
% - Reads the BOLD data from the folders
% - Computes the BOLD phases
% - Calculates the Order Parameter, Mestastability and Synchrony
% - Calculates the instantaneous BOLD synchronization matrix
% - Saves the instantaneous Leading Eigenvectors
%
% Saves the Leading_Eigenvectors and FCD matrices into LEiDA_data.mat
%
% Written by Joana Cabral, May 2016 (joana.cabral@psych.ox.ac.uk)
% Modified by Joana Cabral September 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subjects=9;
Conditions=1:4; % 3=baseline, 4=psilo
TR = 3;

%load data
load psilotc_total.mat tc_aal
[N, Tmax]=size(tc_aal{1,1});

% Create empty variables to save patterns 
Leading_Eig=zeros(Subjects*Tmax,N);
Phase_BOLD=zeros(N,Tmax);
Time_sessions=zeros(2,Subjects*Tmax); % first row = task, 2nd row = subject 

% FILTER SETTINGS             
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

t_all=0;

for task=1:length(Conditions)
    
    for s=1:Subjects
            
        signal = tc_aal{s,Conditions(task)};
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N
            ts=demean(detrend(signal(seed,:)));
            signal_filt =filtfilt(bfilt,afilt,ts);
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
                
        for t=1:Tmax
            t_all=t_all+1;
            %Calculate the Instantaneous FC (BOLD Phase Coherence)
            iFC=zeros(N);
            
            for n=1:N
                for p=1:N
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
            
            % Get the leading eigenvector
            val1=eigs(iFC,1);

            Leading_Eig(t_all,:)=val1;

            Time_sessions(:,t_all)=[task s];
        end
    end
end

save('LEiDA_psilo_data','Metasta','Synchro','Leading_Eig','Time_sessions')