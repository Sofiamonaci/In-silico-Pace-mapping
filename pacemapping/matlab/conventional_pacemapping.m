

function corr = conventional_pacemapping(vt_file,paced_file,flag,out,type_sig,cycle)

%% Usage
% corr = conventional_pacemapping(vt_file,paced_file,flag,out, varargin)
% Compute correlation maps between paces and a reference VT. Signals can be
% either 12-lead ECGs or 8-lead EGMs, however default is 12-lead ECGs, so
% enter 5th argument to compute correlation between 8-lead EGMs
%
% Inputs:
% vt_file:      12-lead ECG or 8-lead EGM VT file/array
% paced_file:   COMBINED pacing 12-lead ECG or 8-lead EGM file/array
%               size: [n_paces*leads, n_time]
% flag:         0: mean of 12/12 or 8/8 correlations
%               1: mean of 10/12 or 6/8 correlations
%               2: Highest 2 leads
% out:          output file (if provided) with .dat extension
% type_sig:     'ECG' or 'EGM'
% cycle:        bcl/cycle of vt_file
%
% Output:
% corr:         correlation map on pacing cloud point
%
%
% Sofia Monaci
% 11/10/21

clc;

fprintf('\n\nCOMPUTING CONVENTIONAL PACE-MAPs FOR VT LOCALISATION...\n\n');

% Loading and Reading files
if isa(vt_file,'char') || isa(vt_file,'string')
    fprintf(' Reading %s ... \n',vt_file);
    vt = dlmread(vt_file,',',0,0);
    vt = vt(:,1:10:end);
else
    vt = vt_file;
end

if isa(paced_file,'char') || isa(paced_file,'string')
    fprintf(' Reading %s ... \n',paced_file);
    pace = dlmread(paced_file,',',0,0);  
else
    pace = paced_file;
end

% Deciding whether to deal with ECGs or EGMs
if contains(type_sig, 'EGM')
    fprintf('Considering 8-lead EGMs ...\n');
    N_leads = 8;
else
    fprintf('Considering 12-lead ECGs ...\n');
    N_leads = 12;
end

% Computing number of pacing locations
N_sites = round(size(pace,1)/N_leads);

% Computing correlations between QRSs
l = 0;
corr = zeros(N_sites,1);

for i = 1:N_sites
    
    fprintf('Pace: %d\n\n',i)
    new_pace = pace(1 + l : N_leads + l,:); 
    corr(i) = QRS_align(vt,new_pace,flag,cycle);
    l = l + N_leads;  
end

if ~isempty(out)
    % Printing out file
    fprintf('Printing out correlation maps for pacing cloudpoints in %s ...\n',out);
    dlmwrite(out,corr*100);
end

end
