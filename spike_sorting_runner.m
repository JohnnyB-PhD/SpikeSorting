% How to perform draft sorting with wave_clus
% (1) Create a new folder with a relevant name.
% (2) Copy all .nev files to process together to the created folder. 
%       The .nev files must be the stripped version 
%       (i.e., processed with "NEV File Information Viewer.exe"). 
% (3) In the space below, add a new text line like the following, 
%           pipeline_for_wave_clus('C:\ephys\186_3_23\188_3_23_LTMFK')
%       where the full path to the created folder (wrapped with single
%       quotes) should appear inside the parentheses. It is possible to 
%       add more than one line for batch processing. Make sure the lines 
%       from past runs are deleted or commented out 
%       (i.e., the line starts with "%" symbol). 
% (4) In MATLAB, open set_parameters.m file in "MATLAB\wave_clus" folder 
%       and verify the sorting parameters.  
% (5) Press Run button in the tool bar of MATLAB. 
%
% When being called, the script does the following.
% (1) Create a new folder named "wave_clus" under the data folder
% (2) Created "*_spikes.mat" files from .nev files in wave_clus folder. 
%       These will be input to the spike sorter.  
% (3) Run the spike sorter. It will create many files including 
%       "times_*.mat" files, which contains the sort results. 
%       This takes tens of minutes to hours, depending on the number of 
%       spikes in the data and sorting parameters. 
% (4) Make new .nev files by combining the sort results and the source 
%       .nev files in the data folder (not wave_clus folder). 
%       The new files have the filenames ends with "rc.nev". 
% 


% pipeline_for_wave_clus('C:\ephys\186_3_23\188_3_23_LTMFK')
pipeline_for_wave_clus('C:\ephys\186_3_23\188_3_23_LTMFK')
