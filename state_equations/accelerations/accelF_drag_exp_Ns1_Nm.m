function [a_drag,varargout] = accelF_drag_exp_Ns1_Nm(state,params,varargin)
% Automatically selects correct exponential drag acceleration
narginchk(2,4)

persistent accelFFun

if ~isempty(accelFFun)
    if nargin == 4
        accelString = "accelF_drag_exp";
        prefix = sprintf("_%ds1",varargin{1});
        suffix = sprintf("_%dm",varargin{2});
    
        accelString = accelString + prefix + suffix;
        accelFFun = str2func(accelString);
    end
else
    if nargin == 4
        accelString = "accelF_drag_exp";
        prefix = sprintf("_%ds1",varargin{1});
        suffix = sprintf("_%dm",varargin{2});
    
        accelString = accelString + prefix + suffix;
        accelFFun = str2func(accelString);
    else
        eid = "Input:missingStateOrModelStateLengths";
        msg = "Please specify both state and model state lengths.";
        error(eid,msg)
    end
end

if nargin == 2
    [a_drag,varargout] = accelFFun(state,params);
end

end