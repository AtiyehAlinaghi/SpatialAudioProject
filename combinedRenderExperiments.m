
% To regenerate the RIR based on the Combined Method



fc = 1000*2.^(-4:4); % central frequencies of the octave bands
% fc = [62.5, 125, 250, 500, 1000, 2000, 4000, 8000,16000];
nLS = 64; % Number of Loudspeakers
[input_stimulus, fs_in] = audioread('speech__HarvardMale.wav');


addpath JSONs\SignalTmix\
addpath RIR
addpath lsir_nondiff_mat_files\

% Load JSON files
cd JSONs\SignalTmix\
jsonFiles = dir('*.json');
cd ..\..
NumOfRooms = length(jsonFiles);

for i = 1:NumOfRooms % Number of rooms

    fileName = jsonFiles(i).name;
    repeat = length('_SignalTmix.json');
    % SignalTmix.json which is repeated for all the json files in our data
    mainPart = fileName(1:(end-repeat));
    % This is for the sake of experiment. If we know the initial delay the
    % BFormat signal is not required
    [rir_gt, fs] = audioread([mainpart '_BFormat.wav']);
    T = round(length(rir_gt)/fs);
    initial_delay = find(abs(rir_gt(:,1))==max(abs(rir_gt(:,1))));
    nondiff_mat = (['lsir_ndiff_' mainPart '.mat']);

    [combined_rir_ls_feed] = combinedRsaoSirrRender(fileName,nondiff_mat, initial_delay, fs,fc,nLS);

    combined_rir_ls_feed = combined_rir_ls_feed.';
    ir_length = size(combined_rir_ls_feed, 1);
    nLS = size(combined_rir_ls_feed, 2);
    output = ...
        matrixConvolver(input_stimulus, reshape(combined_rir_ls_feed, [ir_length, 1, nLS]), size(input_stimulus,1));
    output = 0.99.*output./max(abs(output(:))); % normalise

    FolderLS = ['LS feed ' mainPart];
    cd (FolderLS);

    FolderLS_combined = [mainPart 'combined_LS_feed'];
    status_combined = mkdir(FolderLS_combined);

    cd (FolderLS_Combined);

    for ch =1:nLS
        audiowrite([mainPart num2str(ch) '_combined_LS.wav'], real(combined_rir_ls_feed(:,ch)), fs,'BitsPerSample', 32);
    end

    cd ..

    FolderLS_Speech_combined = [mainPart '_combined_speech_LS_feed'];
    status_combined_speech = mkdir(FolderLS_Speech_combined);

    cd (FolderLS_Speech_combined);


    for ch =1:nLS
        audiowrite([mainPart num2str(ch) '_combined_speech_LS.wav'], real(output(:,ch)), fs,'BitsPerSample', 32);
    end

    cd ..\..
end
%..........................................................................




