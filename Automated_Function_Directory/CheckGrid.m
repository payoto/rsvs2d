function [varargout]=CheckGrid(varargin)
% include_GridCheck
global CheckGrid_Handle
try
nOut=nargout(CheckGrid_Handle);
catch
include_GridCheck
nOut=nargout(CheckGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckGrid_Handle(varargin{:});
end
