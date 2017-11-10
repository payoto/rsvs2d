function [varargout]=MovingAverageCentral(varargin)
% include_Utilities
global MovingAverageCentral_Handle
nOut=nargout(MovingAverageCentral_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MovingAverageCentral_Handle(varargin{:});
end
