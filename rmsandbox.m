function rmsandbox()

fprintf('Removing AstroDyn from path: developer mode\n');

thisFolder = fileparts(mfilename('fullpath'));

rmpath(fullfile(thisFolder, 'tbx'));
rmpath(fullfile(thisFolder, 'tbx', 'astrodyn'));
rmpath(fullfile(thisFolder, 'tbx', 'doc'));

end