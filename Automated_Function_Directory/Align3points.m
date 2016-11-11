function [varargout]=Align3points(varargin)
global Align3points_Handle
nOut=nargout(Align3points_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=Align3points_Handle(varargin{:});
end
