function [varargout]=LengthArea(varargin)
global LengthArea_Handle
nOut=nargout(LengthArea_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LengthArea_Handle(varargin{:});
end
