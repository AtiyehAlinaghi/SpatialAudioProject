% The Script file to read the *.json files containing the RSAO parameters &
% features of each room impulse response (RIR), and use the parameters to
% regenerate the room RIR & the loudspeakers feed to render the signal in
% Unity game engine
%
% Reverberant Spatial Audio Objects (RSAO) algorithm:
% "Object-based reverberation encoding from first-order Ambisonic RIRs"
% Philip Coleman et al. 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Atiyeh Alinaghi, 16/01/2024
%   ECS, University of Southampton
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fc = [62.5, 125, 250, 500, 1000, 2000, 4000, 8000,16000];
fc = 1000*2.^(-4:4); % centre frequencies for octave bands
nLS = 64; % Number of LoudSpeaker
% A sample speech signal
[input_stimulus, fs_in] = audioread('speech__HarvardMale.wav');

%% LOAD JSON FILES

cd JSONs\SignalTmix\
jsonFiles = dir('*.json');
cd ..\..
NumOfRooms = length(jsonFiles);

addpath JSONs\SignalTmix\
addpath RIRs

for i = 2:4 % for Courtyard, DWRC & Kitchen or 1:NumOfRooms
    fileName = jsonFiles(i).name;
    % 'SignalTmix.json' which is repeated for all the json files in our data
    repeat = length('_SignalTmix.json');

    % The parameters file name which is the room name has been chosen to name the
    % folder containing the LoudSpeaker input signals for Unity
    mainPart = fileName(1:(end-repeat));
    % recorded B-Format ground truth (gt) RIR
    [rir_gt, fs] = audioread([mainPart '_BFormat.wav']);
    T = round(length(rir_gt)/fs);

    % To regenerate the early reflections using the parameters in the JSON
    % file. Then VBAP is applied to estimate the LS feeds for rendering part
    ls_setup = load('DTU_ls_dirs_deg.mat');
    [regenEarlyRIR, LS_feed_earlyr] = rsaoEarlyReflectRenderer(fileName,ls_setup,fs,T);

    % To regenerate the Late reverb based on parameters in JSON file. Then
    % decorrelated by random phase for LS feed
    [regenLateReverb, LS_feed_laterev]  = rsaoLateReverbRenderer(fileName,fs,fc,nLS,T);

    LS_feed_rsao_nolag = LS_feed_earlyr+ LS_feed_laterev;
    % As Unity misses the direct sound at t=0 we added some delays
    LS_feed_rsao = [zeros(nLS,10), LS_feed_rsao_nolag];


    rsaoRegenRIR = regenEarlyRIR.' + regenLateReverb;
    L = length(rsaoRegenRIR);
    t = (1:L)/fs;
%     figure; plot(t,rsaoRegenRIR,'r'); xlabel('Time [s]'); ylabel('RIR');
%     title('Regenerated RIR using RSAO approach');
%     xlim([0 0.2]);

    %..........................................................................
    rir_len = size(LS_feed_rsao.', 1);
    nLS = size(LS_feed_rsao.', 2);
    rsaoSpeech_LS_feed = ...
        matrixConvolver(input_stimulus, reshape(LS_feed_rsao.', [rir_len, 1, nLS]), size(input_stimulus,1));
    rsaoSpeech_LS_feed = 0.99.*rsaoSpeech_LS_feed./max(abs(rsaoSpeech_LS_feed(:))); % normalise

    % In order to save the loudspeakers feed as *.wav files
    % The folder name is the room name followed by rsao and loudspeakers feed
    FolderLS = ['loudSpeakers_feed_ ' mainPart '/' mainPart '_rsao_rir_LS_feed'];
    status = mkdir(FolderLS);

    cd (FolderLS);

    for ch =1:nLS
        audiowrite([mainPart '_' num2str(ch) '_rsao_rir_LS.wav'],real(LS_feed_rsao(ch,:)), fs);
    end

    cd ..
    
    FolderLS_Speech_rsao = [mainPart '_rsao_speech_LS_feed'];
    status_rsao = mkdir(FolderLS_Speech_rsao);

    cd (FolderLS_Speech_rsao);

    for ch =1:nLS
        audiowrite([mainPart '_' num2str(ch) '_rsao_speech_LS.wav'], ...
            real(rsaoSpeech_LS_feed(:,ch)), fs,'BitsPerSample', 32);
    end

    cd ..\..
end


