function [varargout]=MakeVideo(varargin)
% include_PostProcessing
global MakeVideo_Handle
nOut=nargout(MakeVideo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeVideo_Handle(varargin{:});
end
