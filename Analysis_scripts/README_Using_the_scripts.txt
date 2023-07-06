%README

% How to use the scripts:
%
% (1) "VKorn_B_TDTR_preprocess_data_for_fit.mlx"
% Pre-process the experimental data with in order to:
%	(a) remove "NaN" values;
%	(b) correct the zero delay time position;
%	(c) correct the phase shift (better use the same value for all the files).
%
% (2) <Your_Parameter_File_Name.m>
% Prepare a parameter file specifying the thermal layers, beam parameters, etc.
% (Samples are available in the "ParameterFiles" folder.
%
% For (3) and (4), place your parameter file to this directory.
%
% (3) "TDTR_Bidirectional_MAIN_SIM.m"
% (ToDo:) check the sensitivity towards the input parameters
%
% (4) "TDTR_Bidirectional_MAIN_FIT.m"
% Fit the experimental data with a model - first by hand (the shape and the position
% of the curve) and then with "auto == 1".
%
% (5) Should be it!

%Include filename and parameters for sample and measurement setup in a
%"Parameter_Example.m" type function. (original version supplied by Cahill has the filename defined in MAIN_FIT)

%TDTR Bidirectional files eat structs. If the measurement data is in
%matrix form, convert it to struct using "TDTR_data_matrix_to_struct". The
%script also has a line of code to shift the time delay zero position to
%correspond with the rise of Vin peak if it has not been successfully
%shifted in the measurement control software.

%If and when the Vout signal of the data is not continuous across tdelay =
%0, there is a phase shift coming from the instrumentation that needs to be
%fixed. This is done by running the data through "AutoSetPhase" function. Files with "_SHIFTED" 
name extension have already been processed.

%After the phase has been fixed TDTR Bidirectional MAIN FIT can be used.