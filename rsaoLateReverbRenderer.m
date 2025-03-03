
function [rirReverb, ls_feed_late, Tmix] = rsaoLateReverbRenderer(metafile,fs,fc,nLS,T)

% This function is based on the paper:   
% P. Coleman, A. Franck, D. Menzies, and P. J. B. Jackson, “Object-based
% reverberation encoding from first-order ambisonic rirs," Journal of the
% Audio Engineering Society, May 2017., to render the late reverberation.

% This function takes the metadata json file contaning the RIR 
% parameters encoded by RSAO method and using those parameters generates 
% the loudspeakrs's feed to render the late reverberation
%..........................................................................
% INPUT ARGUMENTS
%   metafile                   : a *.json file containing the metadata
%                                parameters for the room impulse respose
%                                based on RSAO approach
%   fs                         : sampling frequency
%   fc                         : a vector containing the octave bands
%                                cnetre frequencies
%   nLS                        : number of loud speakers
%   T                          : length of the room impulse response
%..........................................................................
% OUTPUT ARGUMENTS
% rirReverb                    : the regenerated late reverberation before
%                                decorrelation for each loudspeaker 
% ls_feed_late                 : a TxnLS matrix containing the loudspeakers
%                                feeds for late reverbs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Atiyeh Alinaghi, 15/01/2024
%   ECS, Southampton University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath 'jsonlab-master' 'JSONs'

repeat = length('_SignalTmix.json');
mainPart = metafile(1:(end-repeat));

% To decode a JSON file into MATLAB data.
metadata = loadjson(metafile);
disp (metadata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To  load the late reverb parameters from metadata json file
% which is imported as a 1x1 struct

delay = str2double(metadata.room.lreverb.delay);
% "delay" is the time between the direct sound & the 1st reflection

% The "level" is the peak value Pb as in eq.(2) (P. Coleman et al. 2017)
vectLevel = metadata.room.lreverb.level; % An array of chars
subLevel = str2double(split(vectLevel));
% A vector of double values of level for each octave subband

vectRamp = metadata.room.lreverb.attacktime; % An array of chrs
% A vector of double values of the "rising time" from the 1st early reflection
% up to the mixing time for each octave subband

subRamptime = str2double(split(vectRamp));
% Since the attacktime is the time over which the energy rises as a ramp
% and its after the 1st reflection => mixing time is the delay+attacktime
subMixingtime = delay + subRamptime; 

% since the "attacktime" and "delay" are the same over the subbands, 
% we only consider 1
Tmix = subMixingtime(1);
% The decay curve time constants
vectTc = metadata.room.lreverb.decayconst; % An array of chrs
% A vector of double values of time constants for each octave subband
subTimeconst = str2double(split(vectTc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the slope of the rising line from the 1st reflect up to the Tmix
slope = subLevel./subRamptime; 

% fs = 48000; % 48 kHz

% Number of subbands
nBand = length(subTimeconst);

% T = 10; % 10 ms based on the wav file recordings
t = 0:(1/fs):T-(1/fs); % in ms
% The Linearly rising part of lreverb
linerise = zeros(length(t),nBand);
% The Exponentially Decaying Curve: EDC
edc = zeros(length(t),nBand);

for band = 1:nBand
    edc_part = zeros(length(t),1);
    ts  = (delay:(1/fs):subMixingtime(band)-(1/fs)); % Duration in time of the rising line
    % ts  = (delay:(1/fs):Tmix-(1/fs));
    % Duration in time of the rising line in the case where Tmix is the same for all subbands
    
    % The time length for EDC to drop exponentially by 60 dB from Tmix 
    % (In theory)
    % texp  = subMixingtime(band):(1/fs):-(3*log(10)/subTimeconst(band));

    % The time length for EDC to exponentially drop by 120 dB from Tmix 
    % (to be long enough)
    texp  = subMixingtime(band):(1/fs):-(6*log(10)/subTimeconst(band));

    % The time for EDC from Tmix up to the maximum length of RIR
    % texp  = subMixingtime(band):(1/fs):T; % It's too long

    % To compnesate for the round estimate mismatch
    gap = (round(subMixingtime(band)*fs)-round(delay*fs))-length((ts-delay));
    % The Linearly rising part
    linerise(round(delay*fs)+gap:round(subMixingtime(band)*fs)-1,band) = slope(band)*(ts-delay);

    edc_start = round(subMixingtime(band)*fs);
    edc_stop = edc_start + length(texp);
    
    % Time constants values "decayconst" are negative => There is no need
    % for negative multiplication in the exponential decay
    edc_part(edc_start:edc_stop-1) = subLevel(band)*exp(subTimeconst(band)*(texp-subMixingtime(band)));
    
    if length(edc_part) < length(t)
        
        edc_part (edc_stop:end) = edc_part(edc_stop);

    else
        edc_part = edc_part(1:length(t));
    end

    edc(:,band) = edc_part;
    
    clear ts texp edc_part
end

% The envelope of the late reverb
env=linerise+edc;

% Octave bands centre frequencies
% fc = [62.5, 125, 250, 500, 1000, 2000, 4000, 8000,16000];
% fc = 1000*2.^(-4:4);
% nLS = 64;

rirReverb = synthesizeReverb(fs, fc, env, T);

% To match the accumulated energy from all the loudspeakers to the original
% RIR, it is required to normalize the LS feeds. We found these Max_reverb
% heuristically.
switch mainPart

    case 'Kitchen'
        Max_reverb = 0.07;

    case 'Courtyard'
        Max_reverb = 0.1;

    case 'DWRC'
        Max_reverb = 0.06;

    otherwise
        Max_reverb = input('Please Enter the Max_reverb, the normalization factor:');
end
% Max_reverb = max(abs(rirReverb));

ls_feed_late = zeros(nLS,T*fs); % The Loudspeaker feeds

for i_ls =1:nLS
    randomPhase = 1;
    randomPhi = rand(size(rirReverb))*2*pi-pi;
    if randomPhase
        ls_feed_late(i_ls,:) = (abs(rirReverb)./(Max_reverb*nLS)) .* exp(1i*randomPhi);
    end
end

% One way to normalize the late reverb for LS feed so that the recorded one has the
% same amplitude as the original RIR
% sum_feed = repmat(abs(sum(LS_feed_late,1)),[nLS,1]);
% abs_rir = repmat(abs(rirReverb).',[nLS,1]);
% LS_feed_late_norm = (LS_feed_late./sum_feed).*abs_rir;
% Norm = max(abs(sum(LS_feed_late,1)));
% LS_feed_late_norm = (LS_feed_late/Norm);
% for late reverberation after decorrelation using random-phase filters

%..........................................................................
%%% A different way to introduce random phase to the LS feeds
% % All-pass FIR filter with random phase
% N=512; % filter tap length
%
% for ith_ls = 1:nLS
%     a_phase = pi/2; b_phase = -a_phase;
%     Amp_filt = ones(N,1); % Magnitude response is equal to 1 as it is all-pass
%     rnd_phase = (b_phase-a_phase).*rand(N,1)+a_phase; % Random phase
%     h_filt = Amp_filt.*(cos(rnd_phase)+(sqrt(-1)*sin(rnd_phase)));
%     % Frequency domain coefficients
%     hn = ifft(h_filt); % converted to the time domain
%     LS_feed_late(ith_ls,:) = conv(rir_filt, hn, 'same')/(Max_reverb*nLS);
% end
%..........................................................................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function reverb_rir = synthesizeReverb(fs, fc, env, T)

        %  Simulates an exponential decay reverb tail
        %
        % INPUT ARGUMENTS:
        %
        % fs  :    Sample rate
        % fc  :    Center frequencies of reverberation time bands (octave bands)
        % env :    The envelope of reverb energy at each octave band
        % T   :    length of rir

        % OUTPUT ARGUMENTS:
        %
        % reverb_rir:  The reverberant part of RIR, a Tx1 vector 
        
        
        % number of frequency bands
        % nBands = length(fc);

        time = 0:(1/fs):T-(1/fs); % in ms
        rir_len = length(time);

        % Different random noise for each band
        %         rir_init = randn(rir_len, nBands);

        % Same random noise for each octave band
        rir_init = randn(rir_len, 1);
        % To filter the noise into the octave bands 
        mid_out = filterbank(fc, rir_init,fs);
        rir_band = env.*mid_out;
        reverb_rir = sum(rir_band, 2);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rir_band = filterbank(fc, rir_init,fs)

        % For this function we call a class FilterGeneration from the file
        % FilterGeneration.py in python and apply the methods in that class

        nBands = length(fc);
        rir_band = zeros(length(rir_init),nBands);
        BW = [0.3 ones(1,nBands-1)]; % Bandwidth wich is shorter for lowest
        % centre frequency
        for i = 1:nBands

            filterClass=py.FilterGeneration.FilterGeneration(fc(i),BW(i),fs);

            if i == 1

                coef_low = filterClass.lowpassCoefficientsBW;
                b_coef = cell2mat(cell(coef_low.b));
                a_coef = cell2mat(cell(coef_low.a));
                rir_band(:,i) = filtfilt(b_coef, a_coef, rir_init);
                clear a_coef b_coef filterClass

            elseif i == nBand

                coef_high = filterClass.highpassCoefficientsBW;
                b_coef = cell2mat(cell(coef_high.b));
                a_coef = cell2mat(cell(coef_high.a));
                rir_band(:,i) = filtfilt(b_coef, a_coef, rir_init);
                clear a_coef b_coef filterClass
            else

                coef_mid = filterClass.bandpassCoefficientsBW;
                b_coef = cell2mat(cell(coef_mid.b));
                a_coef = cell2mat(cell(coef_mid.a));
                rir_band(:,i) = filtfilt(b_coef, a_coef, rir_init);
                clear a_coef b_coef filterClass
            end

        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end