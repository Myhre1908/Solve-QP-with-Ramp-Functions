function [Q, c, H_ineq, h_ineq, A_eq, b_eq,c0,solution,selectedFile,main_folder] = get_model_data_from_git(model_index)
    %This is the same function as which I have called get_all_data
    %elsewhere, however I improved the folder structure.
    [mat_files_folder,solution_folderpath,c0_folder_path,main_folder] = get_correct_folders();

    [data,selectedFile] = loadMatFileByIndex(model_index,mat_files_folder);

    solution = get_solution(model_index,solution_folderpath,selectedFile);

    c0 = get_c0(model_index,c0_folder_path,selectedFile);

    [Q, c, H_ineq, h_ineq, A_eq, b_eq] = create_matrices(data);
end

function [mat_files_folder,solution_folderpath,c0_folder_path,main_folder] = get_correct_folders()
    filename = 'README.txt';
    fileContent = fileread(filename);
    % Extract all strings inside single quotes (folder paths)
    folderPaths = regexp(fileContent, '"(.*?)"', 'tokens');
    
    % Convert from nested cell to plain cell array
    folderPaths = [folderPaths{:}];
    
    % Display the results
    % for i = 1:length(folderPaths)
    %     fprintf('Path %d: %s\n', i, folderPaths{i});
    % end
    
    main_folder = folderPaths{2};
    mat_files_folder1 = folderPaths{3};
    solution_folder_path1 = folderPaths{4};
    c0_folder_path1 = folderPaths{5};
    
    mat_files_folder = fullfile(main_folder, mat_files_folder1);
    solution_folderpath = fullfile(main_folder, solution_folder_path1);
    c0_folder_path = fullfile(main_folder, c0_folder_path1);
end

function c0 = get_c0(model_index,folderPath,selectedFile)
    filename1 = "stored_c0.mat";
    filename = fullfile(folderPath, filename1);
    data = open(filename).stored_c0;
    mat_file_name = data(model_index,1);
    c0 = double(data{model_index,2});

    if ~strcmp(mat_file_name, selectedFile)
        error("WRONG FILE for c0");
    end
end


function [Q, c, H_ineq, h_ineq, A_eq, b_eq] = create_matrices(data)
    Q = data.Q;  %Quadratic weight
    c = data.c;  %Linear    weight
    A = data.A;  %Equality  matrix
    rl = data.rl; 
    ru = data.ru;
    lb = data.lb;
    ub = data.ub;
    
    %https://github.com/YimingYAN/QP-Test-Problems/blob/master/README.md
    if rl == ru %Ax = b  
        A_eq = A;
        b_eq = rl;
        H_ineq = [speye(length(lb));-speye(length(lb))];
        h_ineq = [ub;-lb];
        % H_ineq = [];
        % h_ineq = [];
        % lb = [];
        % ub = [];
    else %rl \leq Ax \leq ru
        H_ineq = [A;-A;speye(length(lb));-speye(length(lb))];
        h_ineq = [ru;-rl;ub;-lb];
        A_eq = [];
        b_eq = [];
    end
end


function [data,selectedFile] = loadMatFileByIndex(index,folderPath)
    matFiles = dir(fullfile(folderPath, "*.mat"));

    % Check if the folder contains .mat files
    if isempty(matFiles)
        error("No .mat files found in the specified folder.");
    end
    if index < 1 || index > length(matFiles)
        error("Invalid file index. Please choose a number between 1 and %d.", length(matFiles));
    end

    % Load the selected file
    selectedFile = matFiles(index).name;
    data = load(fullfile(folderPath, selectedFile));
    % fprintf("Loaded file: %s\n", selectedFile);
end
function solution = get_solution(index,folderPath, selectedFile)
    % Set the path to the .qp file
    filename1 = '00README.QP';
    filename = fullfile(folderPath, filename1);

    % Open the file
    fid = fopen(filename, 'rt');

    % Skip the header line
    fgetl(fid);

    current_index = 0;
    solution = [];

    while ~feof(fid)
        line = fgetl(fid);
        tokens = regexp(line, '^\s*(\S+)\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+([-+]?\d+(\.\d*)?(e[+-]?\d+)?)\s*$', 'tokens', 'ignorecase');

        if ~isempty(tokens)
            current_index = current_index + 1;

            if current_index == index
                solution = str2double(tokens{1}{2});
                name = upper(tokens{1}{1});
                mat_file_name = [name, '.mat'];
                if ~strcmpi(erase(mat_file_name, '_'), erase(selectedFile, '_'))
                    fclose(fid);
                    error("WRONG FILE for solution");
                end

                break;
            end
        end
    end
    fclose(fid);
    if isempty(solution)
        error('Index out of range.');
    end
end
