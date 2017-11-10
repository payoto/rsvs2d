function [varargout]=ProcDatTemplate(varargin)
% include_Validation
global ProcDatTemplate_Handle
nOut=nargout(ProcDatTemplate_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProcDatTemplate_Handle(varargin{:});
end
