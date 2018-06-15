function [varargout]=GradProcessor(varargin)
% include_Utilities
global GradProcessor_Handle
nOut=nargout(GradProcessor_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GradProcessor_Handle(varargin{:});
end
