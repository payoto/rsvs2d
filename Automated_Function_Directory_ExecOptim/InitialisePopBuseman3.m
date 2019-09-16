function [varargout]=InitialisePopBuseman3(varargin)
% ExecuteOptimisation
global InitialisePopBuseman3_Handle
nOut=nargout(InitialisePopBuseman3_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialisePopBuseman3_Handle(varargin{:});
end
