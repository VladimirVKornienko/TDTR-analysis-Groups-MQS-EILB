%% This script removes all the lines with NaN from the
% processed TDTR files with "*FIN.mat" filenames,
% before running the "TDTR_Bidirectional_MAIN_FIT.m"

% LOAD THE DATA FROM *.MAT-file MANUALLY!
% There should be the <Data> structure.
% It is processed here then, and you can save it back to the file.

%% SOLVED:

Data2 = struct2table(Data);
tmpBefore = length(Data.stagePosition);

Data2 = Data2(~any(ismissing(Data2),2),:);
Data2 = table2struct(Data2,'ToScalar',true);

tmpAfter= length(Data2.stagePosition);
fprintf("Deleted %d entries (%d - %d).\n", round(tmpBefore-tmpAfter),tmpBefore,tmpAfter);

%% best try:
%Data2 = Data;
%myFieldNames = fieldnames(Data2);
%myDim = length(Data.stagePosition);
%Data3 = struct2array(Data2); % testing
%Data4 = array2table(Data3);

%% Redone (WORKS!!!):
%Data2 = Data;
%Data3 = struct2table(Data2);
%Data4 = Data3(~any(ismissing(Data3),2),:);
%Data5 = table2struct(Data4,'ToScalar',true);

%% WRONG TRIES:

% Data2 = struct2array(Data2);

% Data3 = cell2mat(struct2cell(Data2)); % this will convert your structure to matrix


%Data2 = Data2(all(~isnan(Data2),2),:); % for nan - rows
%Data2 = Data2(:,all(~isnan(Data2)));   % for nan - columns

% did not work: >>>
% Data2(isnan([Data2.stagePosition])) = [];   % erase all lines where stagePosition is NaN
% <<<

% Data4 = mat2cell(Data3,ones(myDim,1)*length(myFieldNames));
% Data4 = mat2cell(Data3,[]);

% Data5 = cell2struct(Data4 , myFieldNames, length(myFieldNames)); % this will convert it again to struct

%% END.