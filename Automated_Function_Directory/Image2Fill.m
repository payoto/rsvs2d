function [varargout]=Image2Fill(varargin)
% include_Utilities
global Image2Fill_Handle
nOut=nargout(Image2Fill_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=Image2Fill_Handle(varargin{:});
end
