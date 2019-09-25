function [varargout]=MovingAverageCentral(varargin)
% include_Utilities
global MovingAverageCentral_Handle
try
nOut=nargout(MovingAverageCentral_Handle);
catch
include_Utilities
nOut=nargout(MovingAverageCentral_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingAverageCentral_Handle(varargin{:});
end
