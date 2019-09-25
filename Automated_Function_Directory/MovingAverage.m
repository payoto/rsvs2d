function [varargout]=MovingAverage(varargin)
% include_Utilities
global MovingAverage_Handle
try
nOut=nargout(MovingAverage_Handle);
catch
include_Utilities
nOut=nargout(MovingAverage_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingAverage_Handle(varargin{:});
end
