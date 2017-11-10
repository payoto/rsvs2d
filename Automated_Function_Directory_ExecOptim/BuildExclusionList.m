function [varargout]=BuildExclusionList(varargin)
global BuildExclusionList_Handle
nOut=nargout(BuildExclusionList_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildExclusionList_Handle(varargin{:});
end
