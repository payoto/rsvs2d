function [varargout]=TestDeviationPopStruct(varargin)
% include_Optimisation
global TestDeviationPopStruct_Handle
nOut=nargout(TestDeviationPopStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TestDeviationPopStruct_Handle(varargin{:});
end
