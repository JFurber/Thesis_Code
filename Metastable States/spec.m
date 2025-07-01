function EV_diff = spec(EV)
EV_diff = [];
for i = 1:length(EV)-1
    EV_diff(i) = abs(EV(i)-EV(i+1));
end
end

% THIS CODE CALCULATES THE DIFFERENCE BETWEEN THE EIGENVALUES TO FIND THE SPECTRAL GAP