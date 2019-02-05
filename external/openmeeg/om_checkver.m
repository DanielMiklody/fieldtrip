function version = om_checkver

% Check if OpenMEEG binaries are installed and work. This function
% will return the prefix required to execute the binaries.
%
% Copyright (C) 2010-2017, OpenMEEG developers

% $LastChangedBy: alegra $
% $LastChangedDate: 2010-09-30 11:15:51 +0200 (Thu, 30 Sep 2010) $
% $Revision$

if ispc
    [status, result] = system('om_project_sensors.exe');
   if status>0
        version=22;
    else
        version=24;
    end
    
else
    % the remainder of the code does not apply to windows
    
    location = {
        ''
        '/usr/bin'
        '/usr/local/bin'
        '/usr/local/openmeeg/bin'
        '/opt/bin'
        '/opt/openmeeg/bin'
        };
    
    % start with an empty return value
    prefix = '';
    
    for i=1:numel(location)
        % check whether the binary can be found
        [status, result] = system(sprintf('which %s', fullfile(location{i}, 'om_project_sensors')));
        if status==0
            prefix = location{i};
            % we found it
            break
        end
    end
    
       if status>0
        version=22;
    else
        version=24;
       end        
 
    
end