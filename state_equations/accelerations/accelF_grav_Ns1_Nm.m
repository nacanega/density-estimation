function [a_grav,varargout] = accelF_grav_Ns1_Nm(state,params,varargin)
% Automatically selects correct exponential drag acceleration
narginchk(2,4)

persistent accelFFun

if ~isempty(accelFFun)
    if nargin > 2
        accelString = "accelF_grav";
        suffix = sprintf("_%ds1_Nm",varargin{1});
    
        accelString = accelString + suffix;
        accelFFun = str2func(accelString);
    end
else
    if nargin > 2
        accelString = "accelF_grav";
        suffix = sprintf("_%ds1_Nm",varargin{1});
    
        accelString = accelString + suffix;
        accelFFun = str2func(accelString);
    else
        eid = "Input:missingStateOrModelStateLengths";
        msg = "Please specify both state and model state lengths.";
        error(eid,msg)
    end
end

if nargin == 2
    [a_grav,varargout] = accelFFun(state,params);
end

end