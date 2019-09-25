function [varargout]=GradProcessor(varargin)
% include_Utilities
global GradProcessor_Handle
try
nOut=nargout(GradProcessor_Handle);
catch
include_Utilities
nOut=nargout(GradProcessor_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GradProcessor_Handle(varargin{:});
end
