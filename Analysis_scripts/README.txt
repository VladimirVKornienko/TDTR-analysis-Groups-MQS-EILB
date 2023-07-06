%README

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