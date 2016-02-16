function [varargout]=EdgeFillInformation(varargin)
global EdgeFillInformation_Handle
nOut=nargout(EdgeFillInformation_Handle);
[varargout{1:nOut}]=EdgeFillInformation_Handle(varargin{:});
end
