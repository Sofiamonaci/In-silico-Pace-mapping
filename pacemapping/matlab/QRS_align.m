
function mean_corr = QRS_align(x,y,flag,cycle,varargin)

%% Usage:
% This function allows to align and crop QRSs, used for computing pace-maps
% x:            12/8 x n_times -> first signal (vt or paced)
% y:            12/8 x n_times -> second signal (paced)
% flag:         type of correlation (0: mean of all, 1: mean of highest 10, 2:
%               Leads EGM)
% cycle:        cycle/bcl of x
% varargin:     plotting on
% Output:
% mean_corr:    mean correlation
%
% Sofia Monaci

inv = 0;

% Initialise sequence plots and subplots if plotting on
if size(x,1)==12
    sequence = {'V1','V2','V3','V4','V5','V6','LI','LII','LIII','aVR','aVL','aVF'};
    h =@(x) subplot(4,3,x);
else
    h =@(x) subplot(4,2,x);
    sequence = {'CAN-SVC','CAN-RVcoil','SVC-RVcoil','RVtip-RVring','LVtip-RVtip','LVtip2-RVtip','LVtip3-RVtip','LVtip4-RVtip'};
end

% Loop over each 12-lead/8-lead signal
corr = zeros(1,size(x,1));
for i=1:size(x,1)
    
    % Align signal based on cross-correlation
    [xa,ya,D] = alignsignals(x(i,:),y(i,:));
    
    % Crop aligned signals so they have same size
    if size(xa,2)>size(ya,2)
        xa = xa(1:size(ya,2));
    else
        ya = ya(1:size(xa,2));
    end
    
    if abs(D) > 0
        xa = xa(abs(D):end);
        ya = ya(abs(D):end);
    end
    
    % If signal is at least double the cycle length (meaning that it has more than two QRSs)
    if length(xa) >= 2*cycle
    
        % Find peaks
        [vt_peaks,~] = findpeaks(xa,'MinPeakDistance',cycle);
        
        % Clip QRSs according to intersections
        if length(vt_peaks)>2
            
            [~,j_max] = max(abs(ya));
            [~,~] = max(abs(xa));
            
            % Find first intersection
            for j=j_max:-1:1
                
                if ya(j)*ya(j_max)<=0.006
                    break;
                end
                
            end
            
            ya = ya(j:end);
            xa = xa(j:end);
            
            if (length(xa) - cycle)<=10
                [vt_peaks,~] = findpeaks(xa,'MinPeakDistance',cycle/2);
            else
                [vt_peaks,~] = findpeaks(xa,'MinPeakDistance',cycle);
            end
            
            if length(vt_peaks)> 2
                
                for j=j_max:size(ya,2)
                    
                    if abs(ya(j)*xa(j))<=0.006 || ya(j)*ya(j_max)<=0.009
                        break;
                    end
                    
                end
                
                ya = ya(1:j);
                xa = xa(1:j);
                
            end
            
        end
    end
    
    % Compute correlation between cropped and aligned signals
    R = cov(xa,ya)./(std(xa)*std(ya));
    
    corr(i) = R(1,2);
    
    if nargin > 4
        h(i);
        plot(xa,'-k');
        hold on;
        plot(ya,'--k');
        title({sequence{i},'Corr: ',num2str(corr(i))});
        legend('measured','simulation')
    end
end

if sum(isnan(corr))>0
    
    i = find(isnan(corr));
    inv = length(i);
    
    if length(inv)>2
        error('Too many naan!!');
    end
    
    for j=1:length(i)
        corr = [corr(1:i-1),corr(i+1:end)];
    end
    
end

% Flag 0 returns mean of all correlation, Flag 1 returns mean of best 10
% correlations, Flag 2 - 2 lead for EGM

if flag==0
    mean_corr = mean(corr);
elseif flag==1
    [B,~] = sort(corr,'descend');
    
    if inv==1
        mean_corr = mean(B(1:end-1));
    elseif inv==2
        mean_corr = mean(B);
    else
        mean_corr = mean(B(1:end-2));
    end
    
elseif flag==2
    fprintf('Calculating mean of 2 leads EGMs - NF/FF: %s and %s\n',sequence{2},sequence{4})
    fprintf('Correlation values -> FF: %d\nNF: %d\n',corr(2),corr(4))
    mean_corr = mean(corr([2,4]));
else
    error('Flag must be between 0 and 2, only integerd accepted!)');
    
end

end
