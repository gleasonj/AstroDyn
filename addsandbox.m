function addsandbox()

fprintf('Adding AstroDyn to path: developer mode\n');

thisFolder = fileparts(mfilename('fullpath'));

addpath(fullfile(thisFolder, 'tbx'));
addpath(fullfile(thisFolder, 'tbx', 'astrodyn'));
addpath(fullfile(thisFolder, 'tbx', 'doc'));

end