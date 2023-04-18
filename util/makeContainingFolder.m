function [] = makeContainingFolder(filepath)
%
%
% utility that makes the containing folder for a specified file


[path, file, ext] = fileparts(filepath);

if ~isfolder(path)
    mkdir(path)
end

