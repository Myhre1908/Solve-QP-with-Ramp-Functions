function Give_time_estimate_for_finding_solution(FileName,Main_folder)
    folder_path_from_main_folder = "Data_files/Ramp_solutions";
    try
        Folder_path = fullfile(Main_folder,folder_path_from_main_folder);
        [~, time, ~] = Load_ramp_solution(FileName,Folder_path);
        if time > 30
            Anounce_time(time)
        end
    catch
    end
end


function [Ramp_solution, time, actset] = Load_ramp_solution(FileName,Folder_path)
    

    % Add suffix to filename before loading
    [~, nameOnly, ~] = fileparts(FileName);
    FileNameFinal = nameOnly + "_RampSolution.mat";

    FullPath = fullfile(Folder_path, FileNameFinal);

    % Load the struct from the MAT file
    loadedStruct = load(FullPath, 'RampData');

    % Extract individual fields
    Ramp_solution = loadedStruct.RampData.Ramp_solution;
    time = loadedStruct.RampData.time;
    actset = loadedStruct.RampData.actset;
end

function Anounce_time(seconds)
    disp("Solving this is going to take approximately:");
    hours = floor(seconds / 3600);                % Calculate hours
    minutes = floor(mod(seconds, 3600) / 60);     % Calculate minutes
    seconds = round(mod(seconds, 60));            % Round seconds to nearest whole number
    if hours ~=0
        fprintf('%d hours, %d minutes, %d seconds\n', hours, minutes, seconds);
    else
        fprintf('%d minutes, %d seconds\n', minutes, seconds);
    end
    disp(" ")
end