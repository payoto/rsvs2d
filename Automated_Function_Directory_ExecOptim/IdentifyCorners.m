function [varargout]=IdentifyCorners(varargin)
global IdentifyCorners_Handle
nOut=nargout(IdentifyCorners_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IdentifyCorners_Handle(varargin{:});
end
