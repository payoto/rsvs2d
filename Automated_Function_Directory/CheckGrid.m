function [varargout]=CheckGrid(varargin)
% include_GridCheck
global CheckGrid_Handle
nOut=nargout(CheckGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckGrid_Handle(varargin{:});
end
