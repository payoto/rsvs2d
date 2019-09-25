function [varargout]=MovingAverageLoop(varargin)
% include_Utilities
global MovingAverageLoop_Handle
try
nOut=nargout(MovingAverageLoop_Handle);
catch
include_Utilities
nOut=nargout(MovingAverageLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingAverageLoop_Handle(varargin{:});
end
