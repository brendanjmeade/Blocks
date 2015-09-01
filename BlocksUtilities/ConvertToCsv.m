function ConvertToCsv
% Convert all files with extension in file_type to .csv files
file_types = {'*.sta', '*.segment', '*.block'};
for i = 1:numel(file_types)
    ConvertAll(file_types{i});
end


function ConvertAll(file_type)
% Look in current director for all files of current file_type
dir_data = dir(file_type);
for i = 1:numel(dir_data)
    Sta = ReadStation(dir_data(i).name);
    fprintf(1, 'Converting %s to %s\n', ...
            dir_data(i).name, [dir_data(i).name, '.csv']);
    struct2csv(Sta, [dir_data(i).name, '.csv']);
end


