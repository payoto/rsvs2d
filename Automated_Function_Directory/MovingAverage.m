function [varargout]=MovingAverage(varargin)
global MovingAverage_Handle
nOut=nargout(MovingAverage_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingAverage_Handle(varargin{:});
end
