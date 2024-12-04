function h = nonspherical_alt(r,rhatz,params)
%nonspherical_alt returns the nonspherical altitude given an orbital radius in
% an equauatorial frame

rE = params.rE;
f = params.f;

% Original
%delta = arcsin(rhatz);
%rC = rE*(1-f) / sqrt(1 - (2*f - f^2)*cos(delta)^2)

% Simplified
rC = rE*(1-f) / sqrt(1 - (2*f - f^2)*(1-rhatz^2));
h = r - rC;

end