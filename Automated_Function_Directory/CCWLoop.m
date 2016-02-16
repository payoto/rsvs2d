function [varargout]=CCWLoop(varargin)
global CCWLoop_Handle
nOut=nargout(CCWLoop_Handle);
[varargout{1:nOut}]=CCWLoop_Handle(varargin{:});
end
