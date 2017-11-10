function [varargout]=InitialiseIterationStruct(varargin)
global InitialiseIterationStruct_Handle
nOut=nargout(InitialiseIterationStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseIterationStruct_Handle(varargin{:});
end
