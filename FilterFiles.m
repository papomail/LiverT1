function [ files ] = FilterFiles(folders,filter)
%FilterFiles: selects specific files (based on filename) from multiple folders recursively
myfilter=['*',filter,'*'];

if ischar(folders)
    files = subdir(fullfile(folders,myfilter));
else
    
    NumberLoaded=numel(folders);
    clear files
    files={};
    
    for ii=1:NumberLoaded
        files_1 = subdir(fullfile(folders{ii},myfilter));
        files=cat(2,files,{files_1.name});
    end
    
end