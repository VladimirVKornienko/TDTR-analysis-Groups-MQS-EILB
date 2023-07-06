%% Manual phase-setting subfunction. Same basic functionality as original
%% Adapted from <AutoSetPhase.m>/<PhaseShift> for simpler use from the analysis scripts.
%% Simple subfunction: execute given phase shift on V(out), V(in).
function [Vin_shifted, Vout_shifted ] = VKorn_PhaseShift(phaseShiftDegToApply,Vin_IN, Vout_IN)
    
    radphase=pi/180*phaseShiftDegToApply; % in radians
    
    % Complex representation of phase shift
    VC=(Vin_IN+sqrt(-1)*Vout_IN)*exp(sqrt(-1)*radphase);
    
    % update V(in), V(out) for phase shift
    Vin_shifted  = real(VC);
    Vout_shifted = imag(VC);
    
end