
function [combined_rir_ls_feed] = combinedRsaoSirrRender(fileName,nondiff_mat,initial_delay,fs,fc,nLS)

% This function regenerates the late reverberation ispired by RSAO method &
% combines with the nondiffuse early reflections from SIRR method
%..........................................................................
% INPUT ARGUMENTS
%   fileName                   : a *.json file containing the metadata
%                                parameters for the room impulse respose
%                                based on RSAO approach
%   nondiff_mat                : a *.mat file containing the nondiffuse
%                                part of early reflections based on SIRR
%                                method
%   initial_delay              : Initial delay in the mat file
%   fs                         : sampling frequency
%   fc                         : a vector containing the octave bands
%                                cnetre frequencies
%   nLS                        : number of loud speakers
%..........................................................................
% OUTPUT ARGUMENTS
%         
% combined_rir_ls_feed         : a TxnLS matrix containing the loudspeakers
%                                feeds for the combined method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Atiyeh Alinaghi, 23/01/2024
%   ECS, Southampton University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% To regenerate the Late reverb based on parameters in JSON file. Then
% decorrelated by random phase for LS feed
[~, LS_feed_laterev, Tmix]  = rsaoLateReverbRenderer(fileName,fs,fc,nLS,T);

% In RSAO approach the reverberant energy increases linearly after 1st
% reflection up to the mixing time, Tmix. However, for the combined
% method, this will add extra energy to the early reflections.
% Therefore, we propose a new approach where the reverberant energy
% increases with an envelope of Gaussian curve, with adjustable
% sharpness (variance).
% [regenLateReverb, LS_feed_laterev, Tmix]  = combinedLateReverbRenderer(fileName,fs,fc,nLS,T);
LS_feed_laterev_combined = LS_feed_laterev;

% To load the non-diffuse part of the RIR from SIRR for early
% reflections
lsir_ndiff = load (nondiff_mat);

T_difuse = Tmix*fs + initial_delay;

% To remove the initial delay and the rest of the signal after T_difuse
lsir_ndiff = lsir_ndiff.'; % From the nondiffuse part of SIRR
LS_feed_earlyr_ndiff = [lsir_ndiff(:,initial_delay:T_difuse)...
    zeros(nLS,(T*fs-T_difuse+initial_delay-1))];


% The early reflections are from SIRR
LS_feed_combo_nolag = LS_feed_earlyr_ndiff+ LS_feed_laterev_combined;
% As Unity misses the direct sound at t=0 we added some delays
combined_rir_ls_feed =[zeros(nLS,10),LS_feed_combo_nolag];


