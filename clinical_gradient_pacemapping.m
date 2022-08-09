



function corr_map = clinical_gradient_pacemapping(paced_file,pts_file,out,type_sig,cycle)


%% Usage
% corr = clinical_gradient_pacemapping(paced_file,pts_file,out,varargin)
% Compute correlation maps between paces (gradient pace-mapping) for ishtmus identification. Signals can be
% either 12-lead ECGs or 8-lead EGMs, however default is 12-lead ECGs, so
% enter 5th argument to compute correlation between 8-lead EGMs
%
% Inputs:
% paced_file:   COMBINED pacing 12-lead ECG or 8-lead EGM file/array
%               size: [n_paces*leads, n_time]
% pts_file:     pacing points file/array in cartesian format
% out:          output file (if provided) 
% type_sig:     'ECG' or 'EGM'
% cycle:        bcl/cycle of paced_file
%
% Output:
% corr:         correlation map on pacing cloud point
%
% THRESHOLD IS 20mm, flag = 1 (mean correlation of 10/12 or 6/8)
%
% Sofia Monaci
% 11/10/21

clc;

addpath('/media/sm18/Seagate Backup Plus Drive/PhD/Scripts/ECG/');

fprintf('\n\nCOMPUTING GRADIENT PACE-MAPs FOR VT ISTHMUS LOCALISATION...\n\n');

% Loading and Reading filesout = [out(1:end-4),'.pts'];
if isa(pts_file,'char') || isa(pts_file,'string')
    fprintf(' Reading %s ... \n',pts_file);
    
    if contains(pts_file,'csv')
        pts = dlmread(pts_file,',',0,0); 
    elseif contains(pts_file,'pts')
        pts = dlmread(pts_file,'',1,0);
    else
        error('Pacing points format not compatible! .csv or .pts!');
    end
else
    pts = pts_file;
end

if isa(paced_file,'char') || isa(paced_file,'string')
    fprintf(' Reading %s ... \n',paced_file);
    if contains(paced_file,'.csv')
        pace = dlmread(paced_file,',',0,0);
    elseif contains(paced_file,'.pts')
        pace = dlmread(paced_file,' ',1,0);
    end
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
flag = 1;

% Initialisation variables
count = 1;
corr_map = [];
new_stimuli = [];
m = 0;
thr = 20000; % 20 or 40 mm

for j = 1:N_sites
    
    l = 0;
    pace_old = pace(1 + m : N_leads + m,:);
    pts_old  = pts(j,:);
    
    for i = 1:N_sites
        
        pts_new = pts(i,:);
        pace_new = pace(1 + l: N_leads + l,:);
        
        if i>j && norm(pts_old - pts_new)<thr
            dist = norm(pts_old - pts_new)/1000;
            corr_map(count) = abs(QRS_align(pace_old,pace_new,flag,cycle)*100 - 100)/dist;
            new_stimuli(count,:) = (pts_new + pts_old)/2;
            count = count + 1;
        end
        l = l + N_leads;
    end
    m = m + N_leads;
end

if ~isempty(find(isnan(corr_map), 1)) 
    fprintf('Removing Nan ...\n')
    new_ind = setdiff(1:length(corr_map),find(isnan(corr_map)));
    corr_map = corr_map(new_ind);
    new_stimuli = new_stimuli(new_ind,:);
end

if ~isempty(find(isinf(corr_map), 1)) 
    fprintf('Removing Inf ...\n');
    new_ind = setdiff(1:length(corr_map),find(isinf(corr_map)));
    corr_map = corr_map(new_ind);
    new_stimuli = new_stimuli(new_ind,:);
end

if ~isempty(out)
    
    if ~contains(out,'.dat')
        out = [out,'.dat'];
    end
    
    % Printing out file
    fprintf('Printing out correlation maps for pacing cloudpoints in %s ...\n',out);
    dlmwrite(out,corr_map(:),'delimiter',' ');
    
    out = [out(1:end-4),'.pts'];
    fprintf('Printing out correlation maps for pacing cloudpoints in %s ...\n',out);
    fid = fopen(out,'w');
    fprintf(fid,'%d\n',size(new_stimuli',2));
    fprintf(fid,'%f %f %f\n',transpose(new_stimuli));
    fclose(fid);

end


end

function [final_pts,final_data] = map(data,pts)

count = 1;
final_pts = []; final_data = [];

for j = 1:size(pts,1)
    
    for i=1:size(pts,1)
        
        if j~=i
        final_pts(count,:) = (pts(j,:)+pts(i,:))./2;
        final_data(count) = data(j,i);
        count = count + 1;
        end
        
    end
end

end