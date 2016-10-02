function [varargout]=InitialisePopBuseman(varargin)
global InitialisePopBuseman_Handle
nOut=nargout(InitialisePopBuseman_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialisePopBuseman_Handle(varargin{:});
end
