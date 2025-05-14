function [Ha, Hb, ha, hb, M_HT, selected_indices] = load_initialization_values_for_semi_QP(selectedFile)
    [~, baseName, ext] = fileparts(selectedFile);
    newFileName = [baseName '_initial_values' ext];
    loadedData = load(newFileName, 'myStruct');

    % Extract the fields from the struct
    myStruct = loadedData.myStruct;
    Ha = myStruct.Ha;
    Hb = myStruct.Hb;
    ha = myStruct.ha;
    hb = myStruct.hb;
    M_HT = myStruct.M_HT;
    selected_indices = myStruct.indexes_for_selecting_Ha;
end
