function [varargout]=SortVecColumn(varargin)
% include_Utilities
global SortVecColumn_Handle
try
nOut=nargout(SortVecColumn_Handle);
catch
include_Utilities
nOut=nargout(SortVecColumn_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SortVecColumn_Handle(varargin{:});
end
