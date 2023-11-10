%% Simple subfunction: execute given phase shift on V(out), V(in).
% 10.11.2023.
function [TDTR_Data] = ImportTDTRdataFromFile(filename)
    
    TDTR_Data = struct([]);     % assign empty structure
    
    data = load(filename);      % read data from file
    data = data.Data;

    if isstruct(data)   % preferred option, mostly used
        structData = struct();
        structData.stagePosition = data.stagePosition;%.Data(:,1);
        structData.tdelay = data.tdelay; %.Data(:,2); % in ps
        structData.Vin = data.Vin; %.Data(:,3); %in uV
        structData.Vout = data.Vout; %.Data(:,4); % in uV
        structData.Ratio = data.Ratio; %.Data(:,5);
        structData.Vdet = data.Vdet;
    
        TDTR_Data = structData;
    
    elseif ismatrix(data)   % not tested heavily
        TDTR_Data.stagePosition = data(:,1);
        TDTR_Data.tdelay = data(:,2);
        TDTR_Data.Vin = data(:,3);
        TDTR_Data.Vout = data(:,4);
        TDTR_Data.Ratio = data(:,5);
        TDTR_Data.Vdet = ones(length(data.Data(:,5)),1);
    end
    
end