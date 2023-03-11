function pipeline_for_wave_clus(directory, inputfilename, outputfilename, amp_thr)

% convert a .NEV file to .mat file compatible with wave_clus
%
% Note: this function/script requires the modified versions of openNEV.m 
% and saveNEVSpikes.m, which are named as openNEV_mod.m and 
% saveNEVSpikes_mod.m, respectively. 
%
% [Input arguments] 
% directory: path to the directory caintaining .nev file(s) to process.
%            No backslash ("\") should be added at the end.
%            This input argument cannot be omitted. 
% inputfilename: can be a single filename or blank (i.e., ''). 
%                In the latter case, all .nev files in the directory will 
%                be processed together. 
%                If omitted, it is treated as blank. 
% outputfilename: body of filename for the final outputs. 
%                 If omitted, the name of the directory is used.
% amp_thr: (optional) threshold for eliminating artefacts with very large 
%          amplitudes. Spike events with trough-to-peak amplitude is 
%          (amp_thr * S.D.) will be removed from the dataset. 
% 
% 
% The code can be executed as a script if,
% (1) the filrst line (function definition) is commented out.  
% (2) three variables "workdir", "datafn", and "outputfn" are defined 
%     in the if statement in the beginning. The variables correspond to 
%     "directory", "inputfilename", and "outputfilename", respectively. 
% 

%%
% [1] setting up file information

if ~exist('directory','var')   % when the filrst line is commented out 
                               % and being executed as a script    
clear variables    

% % workdir = 'C:\ephys\Ryota1-4\testR1';
% % datafn = 'R1_103336_Hab1_3_1002.nev';
% % outputfn = '';

% % workdir = 'C:\ephys\186_3_23\wc_test';
% % datafn = '';
% % outputfn = '186_LTM'; 
% % workdir = 'C:\ephys\Ryota1-4\testR4\R4_NLSTM1wc';
% % datafn = '';
% % outputfn = 'R4_NLSTM1wc'; 

% % workdir = 'C:\ephys\Ryota1-4\testR2\R2_NLSTM1wc';
% % datafn = '';
% % outputfn = 'R2_NLSTM1wc'; 

% % workdir = 'C:\ephys\Ryota1-4\testR1\R1_NLSTM1wc';
% % datafn = '';
% % outputfn = 'R1_NLSTM1wc'; 

% % workdir = 'C:\ephys\377\377_0907wc\';
% % datafn = '';
% % outputfn = 'TGStudy_377_9_7wc'; 

% % workdir = 'C:\ephys\Ryota1-4\testR4\R4_NLSTM2wc';
% % datafn = '';
% % outputfn = 'R4_NLSTM2wc'; 

% % workdir = 'C:\ephys\Ryota1-4\testR2\R2_NLSTM2wc';
% % datafn = '';
% % outputfn = 'R2_NLSTM2wc'; 

% % workdir = 'C:\ephys\Ryota1-4\testR1\R1_NLSTM2wc';
% % datafn = '';
% % outputfn = 'R1_NLSTM2wc'; % filename for the output file

% workdir = 'C:\ephys\Ryota1-4\testR4\R4_NLSTM3wc';
% datafn = '';
% outputfn = 'R4_NLSTM3wc'; 

% workdir = 'C:\ephys\Ryota1-4\testR2\R2_NLSTM3wc';
% datafn = '';
% outputfn = 'R2_NLSTM3wc'; 

% workdir = 'C:\ephys\Ryota1-4\testR1\R1_NLSTM4wc';
% datafn = '';
% outputfn = 'R1_NLSTM4wc'; % filename for the output file



else % when being executed as a function
    workdir = directory;
     
    if ~exist('inputfilename', 'var') 
        datafn = '';
    else
        datafn = inputfilename;
    end
    
    if ~exist('outputfilename', 'var') 
        outputfn = '';
    else
        outputfn = outputfilename;
    end
end % if ~exist(directory,'var')
    


% outliers with the trough-to-peak amplitude above amp_thr*S.D. will be discarded. Set 0 to disable.
if ~exist('amp_thr', 'var')
    amp_thr = 6; % this line determines the default value
end

% *** set-up additional parameters ***

sr = 30000; % sampling rate of recording to store in the intermediate files

if ~contains(path, [pwd '\wave_clus'])
    error('Add "wave_clus" folder and its subfolders to MATLAB path.'); 
end

outputdir = [workdir '\wave_clus'];

if ~isempty(datafn)
    file_list = dir([workdir '\' datafn]);
else
    file_list = dir([workdir '\*.nev']);
    if isempty(file_list)
        error('No .nev file in the specified directory.')
    end
    
    if length(file_list) == 1
        disp(['Processing the NEV file "' file_list.name '"...'])
    else
        disp(['Processing ' num2str(length(file_list)) ' NEV files together...'])
    end
end

%% [2] preparation of data files for wave_clus
% reading files

s1 = struct('Electrode', nan, 'TimeStamp', [], 'Waveform', [], 'period', [], 'source_file', []);
% % % all_ch_list = [];

for n = 1:length(file_list)
    
    if isempty(outputfn) && n == 1
        % when outputfn is not specified, use the name of folder that contains the data file(s) 
        outputfn = split(workdir, {'\', '/'});
        outputfn = char(outputfn(end));
        
        % alternative routine: 
        % when outputfn is not specified, take the filename of first data file to be processed 
% %         [~,outputfn] = fileparts(file_list(1).name);

    elseif ~isempty(outputfn)
        [~,outputfn] = fileparts(outputfn);
    end
    
    % requires the modified version of openNEV.m to retrive the all 48 data
    % points of waveforms. The original openNEV.m gets only 46 data points. 
    data=openNEV_mod([workdir '\' file_list(n).name], 'nosave');

    Electrode = data.Data.Spikes.Electrode; 
    ch_list = unique(Electrode); 
% % %     all_ch_list = unique([all_ch_list ch_list]); % currently not used
    
    
    for ch_ct = 1:length(ch_list)
        ch = ch_list(ch_ct);

        flag = (Electrode == ch);

        index = double(data.Data.Spikes.TimeStamp(flag));
        spikes = data.Data.Spikes.Waveform(:,flag)';

        % amplitude filtering
        if amp_thr > 0
            amp = max(spikes,[],2)-min(spikes,[],2);
            stdamp = std(double(amp));
            
            flag2 = amp > amp_thr*stdamp;

            spikes = spikes(~flag2,:);
            index = index(~flag2);
            
        end
        
        s1(n,ch).Electrode = ch;
        s1(n,ch).TimeStamp = index;
        s1(n,ch).Waveform = spikes;
        s1(n,ch).period = [0 data.Data.Spikes.TimeStamp(end)];
        if n > 1
            s1(n,ch).period = s1(n,ch).period + s1(n-1,ch).period(2);
        end
        s1(n,ch).source_file = file_list(n).name;

    end

end

% rearrange the data and generate *_spikes.mat files for wave_clus

mat_file_list = {};
for ch_ct = 1:size(s1,2)
    
    ch = nan;
    for n = 1:size(s1,1)
        if ~isempty(s1(n,ch_ct).Electrode)
            ch = s1(n,ch_ct).Electrode;
            break
        end
    end
    if isnan(ch) % if there is no data associated with the channel number
        continue % skip this channel and go to the next
    end
    
    % assemble output datasets for the current channel
    
    index = [];
    spikes = [];
    period = [];
    source_file = {};
    for n = 1:size(s1,1)
        
        % try to skip when there is no spike in the current channel of
        % current file. Not tested yet. 
        if isempty(s1(n,ch).TimeStamp)
            continue
        end
        
        index = [index (s1(n,ch).TimeStamp + (repmat(double(s1(n,ch).period(1)), [1 length(s1(n,ch).TimeStamp)])))];
        spikes = [spikes; s1(n,ch).Waveform];
        period = [period; s1(n,ch).period];
        source_file = [source_file {s1(n,ch).source_file}];
        
    end
            

    % spike time stamps should be in millisecond in wave_clus
    index = double(index)*1000/sr;
    
    if exist(outputdir, 'dir') == 0
        mkdir(outputdir)
    end
    
% % %     sourceNEV = datafn;
% % %     start_time = 0; % in time point; used in multi-file data
    
    outputpath = [outputdir '\' outputfn '_ch' num2str(ch, '%02.0f') '_spikes'];
    save(outputpath, 'index', 'spikes', 'sr', 'source_file', 'period')
    mat_file_list = [mat_file_list {[outputfn '_ch' num2str(ch, '%02.0f') '_spikes.mat']}];
end

% create the text file for Do_clustering.m
fid = fopen([outputdir '\mat_file_list.txt'], 'w');
for n = 1:length(mat_file_list)

    fprintf(fid, [mat_file_list{n} '\n']);
    
end
fclose(fid);    

% clear(s1) % uncomment this line to remove the full dataset and free up RAM

%% [3] run wave_clus

par = [];

% For running Do_clustering.m with a subset of paramters modified,
% define those paramteres as in the following examples.  
% 
% Default paramteres are defined in "set_paramters.m" file under the 
% "wave_clus" folder. 
%
% [Recommended method]
% If a copy of set_paramters.m is found in the data folder (that contains 
% NEV files and is specified as "directory" argument when called), it will  
% be used as the paraters file rather than the default file in "wave_clus" 
% folder. This method is recommended for running with a custom parameters. 
%
%
% [Obsolete method - not recommended]
% When a parameters is defined here, they will replace the default 
% parameters. 
%
% par.SWCycles = 5000; 
% Increasing this value appears to make the results more stable 
%
% par.min_clus = 50; 
% This value appears to be sensitive to the final number of clusters.
% larger value will result in fewer clusters. 
%

oldFolder = cd(outputdir);

if exist(fullfile(workdir, 'set_parameters.m'), 'file')
    disp('set_paramters.m is found in the datafile folder and will override the default parameter file.')
    copyfile(fullfile(workdir, 'set_parameters.m'));
end

if ~isempty(par)
    Do_clustering('mat_file_list.txt', 'parallel', true, 'par', par)
else
    Do_clustering('mat_file_list.txt', 'parallel', true)
end    
cd(oldFolder);

%% [4] load the sort results 

% read data
result_file_list = dir([outputdir '\times_*.mat']);

s2 = struct('TimeStamp', [], 'Electrode', [], 'Unit', [], 'Waveform', [], 'source_file', {});

for ch_ct = 1:length(result_file_list)
    
    % retrieve channel number
    [~,ch] = fileparts(result_file_list(ch_ct).name);
    ch = replace(ch,'_spikes', ''); % probably unnecessary; just for safety
    ch = ch(end-1:end);
    ch = str2double(ch);
    
    % retrieve some parameters from the input file
    [~,input_file] = fileparts(result_file_list(ch_ct).name);
    input_file = replace(input_file, 'times_', '');
    input_file = dir([outputdir '\' input_file '*.mat']);
    input_data = load([outputdir '\' input_file.name], 'period', 'source_file');
    
    % load the result file for the current channel
    s_temp = load([result_file_list(ch_ct).folder '\' result_file_list(ch_ct).name], 'spikes', 'cluster_class', 'par');
    
    allTimeStamp = round(s_temp.cluster_class(:,2)*s_temp.par.sr/1000);
    allWaveforms = s_temp.spikes;
    
    for n = 1:size(input_data.period,1)
        
        flag = (allTimeStamp > input_data.period(n,1) & allTimeStamp <= input_data.period(n,2));
        
        s2(n,ch).TimeStamp = uint32(allTimeStamp(flag))'-input_data.period(n,1);
        s2(n,ch).Unit = uint8(s_temp.cluster_class(flag,1))';
        s2(n,ch).Electrode = repmat(uint16(ch), [1 length(s2(n,ch).Unit)]);
        s2(n,ch).Waveform = allWaveforms(flag,:)';
        s2(n,ch).source_file = input_data.source_file{n};

    end
end

    
%% [5] rearrange results and reconstruct .NEV files 

for n = 1:size(s2,1)
    
    spikeStruct = struct('TimeStamp', [], 'Electrode', [], 'Unit', [], 'Waveform', []);
    
    for ch_ct = 1:size(s2,2)
        spikeStruct.TimeStamp = [spikeStruct.TimeStamp s2(n, ch_ct).TimeStamp];
        spikeStruct.Unit = [spikeStruct.Unit s2(n, ch_ct).Unit];
        spikeStruct.Electrode = [spikeStruct.Electrode s2(n, ch_ct).Electrode];
        
        spikeStruct.Waveform = [spikeStruct.Waveform s2(n, ch_ct).Waveform];
    end
    
    % rearrange the output array by TimeStamp
    [~,I] = sort(spikeStruct.TimeStamp, 'ascend');
    spikeStruct.TimeStamp = spikeStruct.TimeStamp(I);
    spikeStruct.Electrode = spikeStruct.Electrode(I);
    spikeStruct.Unit = spikeStruct.Unit(I);
    spikeStruct.Waveform = spikeStruct.Waveform(:,I);
    
    % reconstruct NEV file
% %     if n == 1
% %         disp('NOTE: the base .nev file for the output must be "stripped" version of .nev file.')
% %     end
    datafn = s2(n,1).source_file;
    
    NEV = openNEV_mod(fullfile(workdir,datafn), 'nomat', 'nosave');
    
    if ~isempty(NEV.Data.Tracking)
        disp([datafn ': the file appears to be a pre-stripped NEV. Removing the tracking data.'])
        
        NEV.MetaTags.PacketCount = NEV.MetaTags.PacketCount - 2*length(NEV.Data.VideoSync.TimeStamp) - 8;
        NEV.MetaTags.HeaderOffset = NEV.MetaTags.HeaderOffset - 320; 

        NEV = rmfield(NEV, 'VideoSyncInfo');
        NEV = rmfield(NEV, 'ObjTrackInfo');

        NEV.Data.Tracking = [];
        NEV.Data.VideoSync.TimeStamp = [];
        NEV.Data.VideoSync.FileNumber = [];
        NEV.Data.VideoSync.FrameNumber = [];
        NEV.Data.VideoSync.ElapsedTime = [];
        NEV.Data.VideoSync.SourceID = [];
    end
    
    NEV.Data.Spikes.TimeStamp = spikeStruct.TimeStamp;
    NEV.Data.Spikes.Electrode = spikeStruct.Electrode;
    NEV.Data.Spikes.Unit = spikeStruct.Unit;
    NEV.Data.Spikes.Waveform = spikeStruct.Waveform;
    
    
    disp(['Currently saving ' replace(datafn, '.nev', '_rc.nev')]);
    saveNEV_mod(NEV, fullfile(workdir,replace(datafn, '.nev', '_rc.nev')), 'noreport')
% % %     saveNEVSpikes_mod(spikeStruct, replace(datafn, '.nev', '_rc'), workdir, datafn);  
% %     saveNEVSpikes(spikeStruct, replace(datafn, '.nev', '_rc'));  % use this line if saveNEVSpikes_mod.m is not available
    
end

