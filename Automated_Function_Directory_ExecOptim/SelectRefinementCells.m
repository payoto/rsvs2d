function [varargout]=SelectRefinementCells(varargin)
global SelectRefinementCells_Handle
nOut=nargout(SelectRefinementCells_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SelectRefinementCells_Handle(varargin{:});
end
