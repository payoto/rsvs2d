function [varargout]=MovingAverage(varargin)
% include_Utilities
global MovingAverage_Handle
nOut=nargout(MovingAverage_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingAverage_Handle(varargin{:});
end
