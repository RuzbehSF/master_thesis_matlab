
%% Convert Large Excel to MAT Script
% This script reads 'Acceleration_Signals.xlsx' from the 'Real_Data' folder
% and saves a structured .mat file.

close all; clear; clc;

% 1. Setup File Paths
% We use fullfile to ensure it works on both Windows and Mac/Linux
inputFolder = 'Real_Data';
excelFile   = 'Acceleration_Signals.xlsx';
fullPath    = fullfile(inputFolder, excelFile);

outputFile  = fullfile(inputFolder, 'Acceleration_Data.mat');

% 2. Check if file exists
if ~isfile(fullPath)
    error('File not found: %s\nMake sure the script is next to the "Real_Data" folder.', fullPath);
end

fprintf('Detected file: %s\n', fullPath);
fprintf('Analyzing file structure (this may take a moment)...\n');

% 3. Get list of sheets automatically
try
    sheets = sheetnames(fullPath);
    numSheets = length(sheets);
    fprintf('Found %d sheets. Starting import.\n', numSheets);
catch
    error('Could not read sheet names. Ensure you are using a modern version of MATLAB.');
end

% Initialize container structure
AllSignals = struct();

% 4. Loop through sheets and read data
timerVal = tic; % Start timer

for i = 1:numSheets
    currentSheet = sheets(i);
    
    fprintf('Processing Sheet %d/%d: "%s"... ', i, numSheets, currentSheet);
    
    % Read data as a Table
    % 'VariableNamingRule', 'preserve' keeps headers like "sig1_new" exactly as is
    T = readtable(fullPath, 'Sheet', currentSheet, 'VariableNamingRule', 'preserve');
    
    % Create a valid field name for the structure (e.g., "Sheet1")
    % This handles cases where sheet names might have spaces or special chars
    validFieldName = matlab.lang.makeValidName(currentSheet);
    
    % Save table into the structure
    AllSignals.(validFieldName) = T;
    
    fprintf('Done. (%d rows, %d cols)\n', height(T), width(T));
end

% 5. Save to .MAT file
fprintf('\nSaving data to "%s"...\n', outputFile);
fprintf('Note: Using -v7.3 format for large file support.\n');

% Save the structure. '-v7.3' is crucial for big datasets (>2GB)
save(outputFile, 'AllSignals', '-v7.3');

elapsedTime = toc(timerVal);
fprintf('conversion complete in %.2f seconds.\n', elapsedTime);

%% How to use the data later:
fprintf('\n-------------------------------------------------\n');
fprintf('HOW TO LOAD DATA LATER:\n');
fprintf('1. Run: load(''%s'');\n', outputFile);
fprintf('2. Access Sheet 1: data = AllSignals.Sheet1;\n');
fprintf('3. Access specific column: col = AllSignals.Sheet1.sig1_new;\n');
fprintf('-------------------------------------------------\n');