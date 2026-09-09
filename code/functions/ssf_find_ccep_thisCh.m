%% 
% Find ccep

% This function look for an electrode in the file ...task-ccep..events.tsv,
% The output is the number of trials (times that this electrode was
% estimulated)

function [ccepExist, runccepExist] = ssf_find_ccep_thisCh(electrodeName, patientNum, localDataPath)
    %EXAMPLE
    %electrodeName = 'LF2'
    %patientNum = '01';
    
    %To avoid incomplete name. ASK about the second part -
    electrodeName1 = append('-', electrodeName, '	');
    electrodeName2 = append(electrodeName, '-');
    % Open the file for reading in text mode.
    
    directory = fullfile(localDataPath,['sub-', char(string(patientNum))],'ses-ieeg01','ieeg');
    files = dir(fullfile(directory, '*_task-ccep_*events.tsv'));
    ccepExist = 0;
    runccepExist = [];
    for k = 1:length(files)
        filename = [directory '/' files(k).name];
        %fileID = fopen(fullfile(localDataPath, patientNum), '*_electrodes.tsv'), 'rt');
        fileID = fopen(fullfile(filename), 'rt');
        % Read the first line of the file.
        textLine = fgetl(fileID);
        lineCounter = 1;
        
        while ischar(textLine) % While the are characteres, i.e. the final line
            %fprintf('%s\n', textLine);

            % Looking for a string in the line, the result is the number of the
            % column, so if the function doesn't find the string, whereString will be
            % empty and the sum zero
            whereString1 = sum(strfind(textLine, electrodeName1));
            whereString2 = sum(strfind(textLine, electrodeName2));
            
            if whereString2 ~= 0 || whereString1 ~= 0
            % Save the result
                ccepExist = ccepExist + 1; 
                runccepExist = [runccepExist; filename];

            end
            % Read the next line.
            textLine = fgetl(fileID);
            lineCounter = lineCounter + 1;
        end

        % Save the final result
        %sprintf('%d',ccepExist)
        
        % Close the file.
        fclose(fileID);
    end
    runccepExist = unique(string(runccepExist));
end