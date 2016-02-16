function [varargout]=LeftMostCorner(varargin)
global LeftMostCorner_Handle
nOut=nargout(LeftMostCorner_Handle);
[varargout{1:nOut}]=LeftMostCorner_Handle(varargin{:});
end
