function atlas_table = ck_GetAtlasTable(file)

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Code", "Name", "Altname"];
opts.VariableTypes = ["double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Name", "Altname"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Name", "Altname"], "EmptyFieldRule", "auto");

% Import the data
atlas_table = readtable(file, opts);


%% Clear temporary variables
clear opts