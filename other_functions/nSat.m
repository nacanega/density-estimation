function sats = nSat(satScen,constData)
%nSat
% INPUTS:
%
% OUTPUT:
%

nSats = constData.nSats;
satNames = constData.satNames;
names = strings(nSats,1);
for i = 1:nSats
    names(i) = satNames + sprintf("_%d",i);
end

oneVec = ones(nSats,1);

a = constData.semimajor(:);
e = constData.eccent(:);
i = constData.inclin(:);
RAAN = constData.RAAN(:);
ap = constData.argPeri(:);
nu = constData.trueAnom(:);

% Correct sizes of parameters if not specified for each satellite explicitly
if isscalar(a)
    a = a*oneVec*1e3; % Convert to meters and change to correct size.
elseif length(a) ~= nSats
    eid = "Size:semimajorAxisMustBeScalarOrLengthNVector";
    msg = "Semimajoraxis must be a scalar or a vector with length equal" + ...
        " to the number of satellites";
    error(eid,msg);
else
    a = a*1e3; % Convert to meters.
end
if isscalar(e)
    e = e*oneVec;
elseif length(e) ~= nSats
    eid = "Size:eccentricityMustBeScalarOrLengthNVector";
    msg = "Eccentricity must be a scalar or a vector with length equal to" + ...
        " the number of satellites";
    error(eid,msg);
end
if isscalar(i)
    i = i*oneVec;
elseif length(i) ~= nSats
    eid = "Size:inclinationMustBeScalarOrLengthNVector";
    msg = "Inclination must be a scalar or a vector with length equal to" + ...
        " the number of satellites";
    error(eid,msg);
end
if isscalar(RAAN)
    RAAN = RAAN*oneVec;
elseif length(RAAN) ~= nSats
    eid = "Size:RAANMustBeScalarOrLengthNVector";
    msg = "Right ascension of ascending node must be a scalar or a vector" + ...
        " with length equal to the number of satellites";
    error(eid,msg);
end
if isscalar(ap)
    ap = ap*oneVec;
elseif length(ap) ~= nSats
    eid = "Size:argumentPeriapsisMustBeScalarOrLengthNVector";
    msg = "Argument of periapsis must be a scalar or a vector with length" + ...
        " equal to the number of satellites";
    error(eid,msg);
end
if isscalar(nu)
    nu = nu*oneVec;
elseif length(nu) ~= nSats
    eid = "Size:trueAnomalyMustBeScalarOrLengthNVector";
    msg = "True anomaly must be a scalar or a vector with length equal to" + ...
        " the number of satellites";
    error(eid,msg);
end

sats = satellite(satScen,a,e,i,RAAN,ap,nu,Name=names);

end

