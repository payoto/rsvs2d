function [varargout]=BuildSymmetryLists(varargin)
global BuildSymmetryLists_Handle
nOut=nargout(BuildSymmetryLists_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildSymmetryLists_Handle(varargin{:});
end
