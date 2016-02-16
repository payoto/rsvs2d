function [varargout]=EdgeProperties(varargin)
global EdgeProperties_Handle
nOut=nargout(EdgeProperties_Handle);
[varargout{1:nOut}]=EdgeProperties_Handle(varargin{:});
end
