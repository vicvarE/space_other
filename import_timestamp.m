function timestamp=import_timestamp
%% Import data from text file
% Script for importing data from the following text file:


%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["camNum", "frameNum", "sysClock", "buffer"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Setup rules for import
opts.ImportErrorRule = "omitrow";
opts.MissingRule = "omitrow";

% Import the data
timestamp = readtable("timestamp.dat", opts);


%% Clear temporary variables
clear opts
end