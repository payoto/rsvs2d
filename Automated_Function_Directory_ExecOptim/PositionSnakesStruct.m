function [varargout]=PositionSnakesStruct(varargin)
global PositionSnakesStruct_Handle
nOut=nargout(PositionSnakesStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PositionSnakesStruct_Handle(varargin{:});
end
