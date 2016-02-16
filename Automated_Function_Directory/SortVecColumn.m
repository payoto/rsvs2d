function [varargout]=SortVecColumn(varargin)
global SortVecColumn_Handle
nOut=nargout(SortVecColumn_Handle);
[varargout{1:nOut}]=SortVecColumn_Handle(varargin{:});
end
