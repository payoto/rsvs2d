function [varargout]=MovingAverageLoop(varargin)
% include_Utilities
global MovingAverageLoop_Handle
nOut=nargout(MovingAverageLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingAverageLoop_Handle(varargin{:});
end
