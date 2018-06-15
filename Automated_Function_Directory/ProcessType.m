function [varargout]=ProcessType(varargin)
% include_Utilities
global ProcessType_Handle
nOut=nargout(ProcessType_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProcessType_Handle(varargin{:});
end
