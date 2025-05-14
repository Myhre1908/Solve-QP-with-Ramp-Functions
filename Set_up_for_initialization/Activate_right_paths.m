function Activate_right_paths(Main_folder,solver)
    solve_ramp_functions_folder = fullfile(Main_folder,"Solve_ramp_equations");
    addpath(solve_ramp_functions_folder);
    initialize_ramp_functions_folder = fullfile(Main_folder,"Set_up_ramp_functions_solvers");
    addpath(initialize_ramp_functions_folder);

    folderPaths = find_right_folders(Main_folder);
    if solver == "qpramp"
        addpath(folderPaths(1));
    elseif solver == "SemiQP"
        addpath(folderPaths(2));
        addpath(folderPaths(3));
        addpath(folderPaths(4));
        
    elseif solver == "SemiQP_ineq"
        addpath(folderPaths(4));
        addpath(folderPaths(5));
        addpath(folderPaths(6));
    elseif solver == "new_QP"
        addpath(folderPaths(4));
    end
end

function folderPaths = find_right_folders(Main_folder)
    Path_to_File_paths = fullfile(Main_folder,"Set_up_for_initialization");
    filename = 'File_paths.txt';
    fullFilePath = fullfile(Path_to_File_paths, filename);
    fileContent = fileread(fullFilePath);

    folderPaths = regexp(fileContent, '"(.*?)"', 'tokens');
    folderPaths = [folderPaths{:}];
    Full_folder_paths = cellfun(@(f) fullfile(Main_folder, f), folderPaths, 'UniformOutput', false);

    % Display the results
    % for i = 1:length(Full_folder_paths)
    %     fprintf('Path %d: %s\n', i, Full_folder_paths{i});
    % end
    folderPaths = string(Full_folder_paths);
end

