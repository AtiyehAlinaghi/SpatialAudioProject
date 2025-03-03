
function [rirReflections, ls_feed_early] = rsaoEarlyReflectRenderer(metafile,ls_setup,fs,T)

% This function is based on the paper:   
% P. Coleman, A. Franck, D. Menzies, and P. J. B. Jackson, “Object-based
% reverberation encoding from first-order ambisonic rirs," Journal of the
% Audio Engineering Society, May 2017., to render the early refelections.

% This function takes the metadata json file contaning the RIR 
% parameters encoded by RSAO method 
% Then recovers the RIRs and then uses VBAP renderer to calculate the 
% loudspeakers' input to regenerate the room impulse response 
%..........................................................................
% INPUT ARGUMENTS
%  metafile                : a json file containing the RSAO RIR parameters
%  fs                      : sampling frequency
%  ls_setup                : the position of the loudspeakers for the 
%                            rendering in the polar coordinate 
%  T                       : the time length of the RIR in seconds
%..........................................................................
% OUTPUT ARGUMENTS
% rirReflections           : regenerated early reflections in RIR 
% ls_feed_early            : an N_ls X T matrix of the loudspeaker inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Atiyeh Alinaghi, 15/01/2024
%   ECS, University of Southampton
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath 'jsonlab-master' 'JSONs' 'Vector-Base-Amplitude-Panning-master'

metadata = loadjson(metafile);
disp (metadata)

% The azimuth, elevation & radius of the direct sound
% radius = str2double(metadata.position.radius); % It is set to 1
direct_az = str2double(metadata.position.az);
direct_el = str2double(metadata.position.el);

% The length of regenerated RIR
% T = 10; % 10 sec
% fs = 48000; % fs = 48 kHz
L = T*fs; % number of samples in recovered RIR
% The number of early reflection peaks 
num_early = size(metadata.room.ereflect,2);
% The number of 2nd order IIR filters to model each reflection spectrum
num_filter = size(metadata.room.ereflect(1).biquadsos,2);


% rir = [1, zeros(1,L-1)]; 
direct = [1, zeros(1,L-1)];
% mid = [1, zeros(1,L-1)];
reflects = zeros(num_early, L);             
each_azim = zeros(1,num_early); each_elev = zeros(1,num_early);

for i = 1:num_early
    % mid = [1, zeros(1,L-1)];
    each_reflect = metadata.room.ereflect (i);
    each_level = str2double(each_reflect.level);
    each_delay = str2double(each_reflect.delay);
    each_azim (i) = str2double(each_reflect.position.az);
    each_elev (i) = str2double(each_reflect.position.el);

    % Each peak is delayed & scaled 
    mid = each_level * ...
        [zeros(1, round(each_delay*fs)-1), 1, zeros(1,L-round(each_delay*fs))];
    y = mid;
    clear mid
    
    % designing an IIR filter using a series of biquad (2nd order) filters
    % The reason to use a series of biquad filters is that 
    % they can be implemented without stability issues  
    num_coefs = size(fieldnames(metadata.room.ereflect(1).biquadsos(1)),1);
    n_order = (num_coefs/2)-1;
    for j = 1:num_filter        
        coef_b = ones(1,n_order+1); coef_a = ones(1,n_order+1);        
        for r = 0:n_order
        coef_b(r+1) = str2double(each_reflect.biquadsos(j).("b"+ num2str(r)));
        coef_a(r+1) = str2double(each_reflect.biquadsos(j).("a"+ num2str(r)));
        end
        y = filtfilt(coef_b, coef_a, y); % To have zero-phase filter
    end
    reflects (i,:) = y;
    clear mid y
end

rirReflections = direct + sum(reflects); 

% The direct sound & early reflections each in a row
direc_reflec = [direct ; reflects];


% Loudspeaker setup for vector based amplitude panning (VBAP) rendering

% ls_dirs_deg = load('DTU_ls_dirs_deg.mat'); 
ls_dirs_deg = ls_setup;
ls_dirs = ls_dirs_deg.ls_dirs_deg; 
ls_groups = findLsTriplets(ls_dirs);
ls_invMtx = invertLsMtx(ls_dirs,ls_groups);
src_dirs = [direct_az direct_el ; each_azim' each_elev'];
gain3D = vbap(src_dirs,ls_groups,ls_invMtx); 
% The gains for each reflection at each loudspeaker 

ls_feed_early = gain3D' * direc_reflec;
% The loudspeakers' feed
