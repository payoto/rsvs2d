function [varargout]=SortVecColumn(varargin)
global SortVecColumn_Handle
nOut=nargout(SortVecColumn_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SortVecColumn_Handle(varargin{:});
end
