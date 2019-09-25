function [varargout]=ProcDatTemplate(varargin)
% include_Validation
global ProcDatTemplate_Handle
try
nOut=nargout(ProcDatTemplate_Handle);
catch
include_Validation
nOut=nargout(ProcDatTemplate_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProcDatTemplate_Handle(varargin{:});
end
