function [Ha, Hb, ha, hb, Ha_inv, selected_indices] = load_initialization_values_for_semi_QP_ineq(selectedFile)
    % Extract the base name and extension from the selected file
    [~, baseName, ext] = fileparts(selectedFile);
    newFileName = [baseName '_initial_values' ext];
    loadedData = load(newFileName, 'myStruct');
    myStruct = loadedData.myStruct;

    % Unpack the fields
    Ha = myStruct.Ha;
    Hb = myStruct.Hb;
    ha = myStruct.ha;
    hb = myStruct.hb;
    Ha_inv = myStruct.Ha_inv;
    selected_indices = myStruct.indexes_for_selecting_Ha;
end
