function [varargout]=PerformIteration(varargin)
global PerformIteration_Handle
nOut=nargout(PerformIteration_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PerformIteration_Handle(varargin{:});
end
