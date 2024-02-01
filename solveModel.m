%  Time points with experimental data
lineageCloneSizeDays = [1,2,3,4,5,6,34,92,185,367];
lineageCloneSizeWeeks = lineageCloneSizeDays/7;

% load experiment data
loadExperimentData

% max clone size
max_clone_size = 300;

% set parameter sweep values
clear r lam rho
r.min = 0.05;
r.skip = 0.05;
r.max = 0.15;
lam.min = 0.9;
lam.skip = 0.2;
lam.max = 1.3;
rho.min = 0.05;
rho.skip = 0.05;
rho.max = 0.15;

% functions for converting idx to paramter value
r_i2v = @(r,idx) r.min+(idx-1)*r.skip;
lam_i2v = @(lam,idx) lam.min+(idx-1)*lam.skip;
rho_i2v = @(rho,idx) rho.min+(idx-1)*rho.skip;
gam_i2v = @(rho,rho_idx,lam,lam_idx) (lam_i2v(lam,lam_idx)*rho_i2v(rho,rho_idx)/(1-rho_i2v(rho,rho_idx)));

% calculate all clone size probabilities unless already stored
disp('========================================')
disp('Calculating all clone size probabilities')
disp('========================================')
file_name = ['all_clone_size_probs_1m_3m_6m_12m_r_' num2str(r.min) '_' num2str(r.skip) '_' num2str(r.max) '_lam_' num2str(lam.min) '_' num2str(lam.skip) '_' num2str(lam.max) '_rho_' num2str(rho.min) '_' num2str(rho.skip) '_' num2str(rho.max) '.mat'];
if exist(file_name, 'file') == 2
    load(file_name);
 else
   all_clone_size_probs = all_clone_size_probabilities(lineageCloneSizeWeeks(7:10),max_clone_size,r,lam,rho);
save(file_name,'all_clone_size_probs')
end
disp('=============================================')
disp('Done calculating all clone size probabilities')
disp('=============================================')


% Troy data
disp('==========================================================')
disp('Troy experimental data: Calculating parameter probabilties ')
disp('==========================================================')
Troy_all_param_prob_score = model_parameter_fit(all_clone_size_probs,{Td{7:10}},r,lam,rho);
disp('===============================================================')
disp('Troy experimental data: Done calculating parameter probabilties')
disp('===============================================================')

disp('=====================================================')
disp('Troy experimental data: Finding best model parameters')
disp('=====================================================')
clear Troy_bestfit
for jj = 1:4
    Troy_bestfit{jj} = best_parameter_fit(squeeze(all_param_probs(:,:,:,jj)),r,lam,rho);
end
Troy_bestfit_overall = best_parameter_fit_overall(Troy_all_param_prob_score,r,lam,rho);
disp('==========================================================')
disp('Troy experimental data: Done finding best model parameters')
disp('==========================================================')


  
% Sox2 data

disp('==========================================================')
disp('Sox2 experimental data: Calculating parameter probabilties ')
disp('==========================================================')
Sox_all_param_prob_score = model_parameter_fit(all_clone_size_probs(:,:,:,2),{TdSox},r,lam,rho);
disp('===============================================================')
disp('Sox2 experimental data: Done calculating parameter probabilties')
disp('===============================================================')

disp('=====================================================')
disp('Sox2 experimental data: Finding best model parameters')
disp('=====================================================')
clear Sox_bestfit
Sox_bestfit = best_parameter_fit(squeeze(all_param_probsSox),r,lam,rho);
disp('==========================================================')
disp('Sox2 experimental data: Done finding best model parameters')
disp('==========================================================')


% Explain results
disp('  ')
disp('===================================')
disp('============= OUTPUTS =============')
disp('===================================')
disp('  ')
disp('============= CLONE SIZE PROBABILITIES =============')
disp('all_clone_size_probs : contains clone size probabilities for clone sizes 1 to 300 for different r, lam, rho, and time points')
disp('                     : 1st index: rr corresponds to r-value r.min+(rr-1)*r.skip')
disp('                     : 2nd index: ll corresponds to lambda-value lam.min+(ll-1)*lam.skip')
disp('                     : 3nd index: gg corresponds to rho-value rho.min+(gg-1)*rho.skip')
disp('                     : 4th index: tt=1 corresponds to 1 month, tt=2 to 3 months, tt=3 to 6 months, tt=4 to 12 months')

disp('  ')
disp('============= TROY RESULTS =============')
disp('Troy_all_param_prob_score : contains the model parameter fit score for different r, lam, rho, and time points')
disp('                          : smaller value is better fit (and can be translated to a probabilty)')
disp('                          : 1st index: rr corresponds to r-value r.min+(rr-1)*r.skip')
disp('                          : 2nd index: ll corresponds to lambda-value lam.min+(ll-1)*lam.skip')
disp('                          : 3nd index: gg corresponds to rho-value rho.min+(gg-1)*rho.skip')
disp('                          : 4th index: tt=1 corresponds to 1 month, tt=2 to 3 months, tt=3 to 6 months, tt=4 to 12 months')
disp('Troy_bestfit : contains the best model parameters and parameter fit score for different time points')
disp('                          : 1st index: tt=1 corresponds to 1 month, tt=2 to 3 months, tt=3 to 6 months, tt=4 to 12 months')
disp('Troy_bestfit_overall : contains the best model parameters and parameter fit score considering all time points')

disp('  ')
disp('============= SOX2 RESULTS =============')
disp('Sox_all_param_prob_score : contains the model parameter fit score for time point 3 months and different r, lam, and rho')
disp('                         : smaller value is better fit (and can be translated to a probabilty)')
disp('                         : 1st index: rr corresponds to r-value r.min+(rr-1)*r.skip')
disp('                         : 2nd index: ll corresponds to lambda-value lam.min+(ll-1)*lam.skip')
disp('                         : 3nd index: gg corresponds to rho-value rho.min+(gg-1)*rho.skip')
disp('Sox_bestfit : contains the best model parameters and parameter fit score for time point 3 months')



%%%%%%%%%%%% FUNCTION DEFINITIONS %%%%%%%%%%%%

function max_prob = best_parameter_fit_overall(all_probs,r,lam,rho)

max_prob.val = -inf;

for rIdx = 1:size(all_probs,1)
    rVal = r.min+(rIdx-1)*r.skip;
    for lamIdx = 1:size(all_probs,2)
        lamVal = lam.min+(lamIdx-1)*lam.skip;
        for gamIdx = 1:size(all_probs,3)
            rhoVal = rho.min+(gamIdx-1)*rho.skip;
            gamVal = rhoVal*lamVal/(1-rhoVal);
            prob = sum(all_probs(rIdx,lamIdx,gamIdx,:));

            if prob > max_prob.val %&& rhoVal >= rho.min && rhoVal <= rho.max
                 max_prob.val = prob;
                 max_prob.rIdx = rIdx;
                 max_prob.lamIdx = lamIdx;
                 max_prob.gamIdx = gamIdx;
                 max_prob.r = rVal;
                 max_prob.lam = lamVal;
                 max_prob.gam = gamVal;
                 max_prob.rho = rhoVal;
             end
        end
    end
end
end

function max_prob = best_parameter_fit(all_param_probs_at_t,r,lam,rho)

max_prob.val = -inf;

for rIdx = 1:size(all_param_probs_at_t,1)
    rVal = r.min+(rIdx-1)*r.skip;
    for lamIdx = 1:size(all_param_probs_at_t,2)
        lamVal = lam.min+(lamIdx-1)*lam.skip;
        for gamIdx = 1:size(all_param_probs_at_t,3)
            rhoVal = rho.min+(gamIdx-1)*rho.skip;
            gamVal = rhoVal*lamVal/(1-rhoVal);
            prob = all_param_probs_at_t(rIdx,lamIdx,gamIdx);
            if prob > max_prob.val %&& rhoVal >= rho.min && rhoVal <= rho.max
                 max_prob.val = prob;
                 max_prob.rIdx = rIdx;
                 max_prob.lamIdx = lamIdx;
                 max_prob.gamIdx = gamIdx;
                 max_prob.r = rVal;
                 max_prob.lam = lamVal;
                 max_prob.gam = gamVal;
                 max_prob.rho = rhoVal;
             end
        end
    end
end
end

function all_param_probs = model_parameter_fit(all_clone_size_probs,Td,r,lam,rho)

num_r = floor((r.max - r.min) / r.skip + 1);
num_lam = floor((lam.max - lam.min) / lam.skip + 1);
num_rho = floor((rho.max - rho.min) / rho.skip + 1);  % Since gamVal is derived from rhoVal

all_clone_size_p = cell(num_r,num_lam,num_rho,length(Td));
all_clone_size_p(:,:,:,:) = all_clone_size_probs;
all_clone_size_probs = all_clone_size_p;


% 3D array to store all 'prob' values
all_param_probs = zeros(num_r, num_lam, num_rho, length(Td));

for rIdx = 1:num_r
    rVal = r.min + (rIdx-1)*r.skip;
    
    % Temporary 2D array to store 'prob' values for this rVal
    temp_param_probs_for_r = zeros(num_lam, num_rho);
    
    for lamIdx = 1:num_lam
        lamVal = lam.min + (lamIdx-1)*lam.skip;
        
        for rhoIdx = 1:num_rho
            rhoVal = rho.min + (rhoIdx-1)*rho.skip;
            gamVal = rhoVal*lamVal/(1-rhoVal);
            
            for tIdx = 1:length(Td)
                temp_param_probs_for_r(lamIdx, rhoIdx, tIdx) = parameter_logprobabilities_nonnomalized(all_clone_size_probs(rIdx,lamIdx,rhoIdx,tIdx), Td{tIdx}.sum);
            end
        end
    end
    all_param_probs(rIdx, :, :, :) = temp_param_probs_for_r;
end
end

function all_clone_size_probs = all_clone_size_probabilities(t,max_clone_size,r,lam,rho)

N = max_clone_size;
num_r = floor((r.max - r.min) / r.skip + 1);
num_lam = floor((lam.max - lam.min) / lam.skip + 1);
num_rho = floor((rho.max - rho.min) / rho.skip + 1);  

% 3D array to store all 'prob' values
temp_clone_size_probs = zeros(num_r, num_lam, num_rho, length(t), N);

for rIdx = 1:num_r
    rVal = r.min + (rIdx-1)*r.skip;
    fprintf(2,'r:%d\n',rVal)
    
    % Temporary 2D array to store 'prob' values for this rVal
    %temp_param_probs_for_r = zeros(num_lam, num_rho);
    % Temporary storage for clone size probabilities inside parfor
    temp_clone_size_probs_for_r = zeros(num_lam, num_rho, length(t), N);
    
    for lamIdx = 1:num_lam
        lamVal = lam.min + (lamIdx-1)*lam.skip;
        fprintf(2,'lambda:%d\n',lamVal)
        
        for rhoIdx = 1:num_rho
            rhoVal = rho.min + (rhoIdx-1)*rho.skip;
            gamVal = rhoVal*lamVal/(1-rhoVal);
            
            [clone_size_prob,~] = clone_size_probabilities(t,lamVal,gamVal,rVal,max_clone_size);
            for tIdx = 1:length(t)
                temp_clone_size_probs_for_r(lamIdx,rhoIdx,tIdx,:) = clone_size_prob{tIdx};
            end
        end
    end
    temp_clone_size_probs(rIdx,:,:,:,:) = temp_clone_size_probs_for_r;
    %all_param_probs(rIdx, :, :) = temp_param_probs_for_r;
end

% Transfer the temporary clone size probabilities to the final cell array
all_clone_size_probs = cell(num_r, num_lam, num_rho, length(t));
for rIdx = 1:num_r
    for lamIdx = 1:num_lam
        for rhoIdx = 1:num_rho
            for tIdx = 1:length(t)
                all_clone_size_probs{rIdx, lamIdx, rhoIdx, tIdx} = squeeze(temp_clone_size_probs(rIdx, lamIdx, rhoIdx, tIdx, :));
            end
        end
    end
end
end

function log_prob = parameter_logprobabilities_nonnomalized(p,f)
p = p{1};
p = p(length(f));
log_prob_values = -f .* log(p);
log_prob = sum(log_prob_values);

end

function [prob,prob_no_correction] = clone_size_probabilities(t,lambda,gamma,r,max_clone_size)

rho = gamma/(lambda+gamma);
tmp = fft_on_boundary(t,lambda,gamma,r);
for idx = 1:length(t)
    fft_values{idx} = real(tmp{idx});    
    prob{idx} =  fft_values{idx}(2:end)/sum(fft_values{idx}(2:end));
    %log_prob = log(prob{idx});
    prob_no_correction{idx} = prob{idx}(1:max_clone_size)';
    [cons_diff,idx_u] = find_linear_exp_trend(prob{idx}(1:max_clone_size),1e-4);
    prob{idx}(idx_u+1:end) = exp(log(prob{idx}(idx_u))+cons_diff*[1:length(prob{idx})-idx_u]);
    prob{idx} = prob{idx}(1:max_clone_size)';
end
end

function [trend,idx_u] = find_linear_exp_trend(v,tol)

trend = NaN;
v = diff(log(v));

% Compute absolute differences
diffs = abs(diff(v));

while isnan(trend)
    
    % Find regions where the difference is less than tolerance
    regions = [0, diffs < tol, 0];
    
    % Compute differences to identify the start and end of regions
    d = diff(regions);
    
    % Identify starts and ends of each region
    starts = find(d == 1);
    ends = find(d == -1);
    
    % Find the longest region
    [lengths, idx] = max(ends - starts);
    
    % Indices of the longest region in the original vector
    longest_region_indices = starts(idx):ends(idx);
    
    % Values of the longest region
    longest_region_values = v(longest_region_indices);
    
    trend = mean(longest_region_values);
    idx_u = ends(idx);
    tol = tol*10;
end
end


function fft_values = fft_on_boundary(t,lambda,gamma,r)
d = 12;
N = 2^d;

%x_values = linspace(0, 2*pi, N);
x_values = (1:N)*2*pi/N;

for idx = 1:length(t)
    values_on_circle{idx} = arrayfun(@(x) H(exp(i*x), t(idx), lambda,gamma,r)*exp(-2*i*pi/N), x_values);
    % Perform FFT on the resulting vector
    H_fft{idx} = fft(values_on_circle{idx});
    fft_values{idx} = H_fft{idx}/N;
end

end

function Hval = H(x,t,lambda,gamma,r)
Hval = G(x,x,t,lambda,gamma,r);
end

function Gval = G(x,y,t,lambda,gamma,r)
% Initial conditions [Fx0, Fy0]
y0 = [x; y];

% Time span
tspan = [0 t];

% Solve ODE
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
[t, Y] = ode45(@(t, z) coupledODE(t, z, lambda, gamma, r), tspan, y0, options);
%[t, Y] = ode15s(@(t, z) coupledODE(t, z, lambda, gamma, r), tspan, y0, options);

Gval = Y(end,1);

end

function dydt = coupledODE(t, y, lambda, gamma, r)
% Extract variables from y
Fx = y(1);
Fy = y(2);

% Define f_x and f_y
fx = r*Fx^2 + (1-2*r)*Fx*Fy + r*Fy^2;
fy = 1;

% Define differential equations
dFx = lambda * (fx - Fx); % Here, I'm assuming 'x' is equivalent to 't'. Please adjust if different.
dFy = gamma * (fy - Fy); % Here, I'm assuming 'y' is equivalent to 't'. Please adjust if different.

dydt = [dFx; dFy];
end
