function [varargout]=InitialisePopBuseman2(varargin)
% ExecuteOptimisation
global InitialisePopBuseman2_Handle
nOut=nargout(InitialisePopBuseman2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialisePopBuseman2_Handle(varargin{:});
end
