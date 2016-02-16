function [varargout]=OrderBlockEdges(varargin)
global OrderBlockEdges_Handle
nOut=nargout(OrderBlockEdges_Handle);
[varargout{1:nOut}]=OrderBlockEdges_Handle(varargin{:});
end
