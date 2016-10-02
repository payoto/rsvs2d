function [varargout]=InitialiseOptimisation(varargin)
global InitialiseOptimisation_Handle
nOut=nargout(InitialiseOptimisation_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseOptimisation_Handle(varargin{:});
end
