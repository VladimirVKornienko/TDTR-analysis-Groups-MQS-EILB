%% Remove every row from the structure <Data> that
%% contains at least one NaN-value
% 10.11.2023.
% <verboseFlag>     : to print the output or not.
function [TDTR_Data_Out] = RemoveNaNfieldsFromTheStruct(TDTR_Data_In, verboseFlag)
    
    TDTR_Data_Out = struct([]);

    %%% Apr. 19 >>> %%%
    % Remove 'NaN' values: %
    tmpBefore = length(TDTR_Data_In.stagePosition);

    TDTR_Data_Out = struct2table(TDTR_Data_In);
    TDTR_Data_Out = TDTR_Data_Out(~any(ismissing(TDTR_Data_Out),2),:);
    TDTR_Data_Out = table2struct(TDTR_Data_Out,'ToScalar',true);
    
    tmpNumDelItems = round(tmpBefore-length(TDTR_Data_Out.stagePosition));
    
    if (verboseFlag)
        if (tmpNumDelItems == 1)
            fprintf("Deleted 1 entry.")
        else
            fprintf("Deleted %d entry(ies).", length(TDTR_Data_Out.stagePosition))
        end
    end
    %%% <<< %%%
    
end