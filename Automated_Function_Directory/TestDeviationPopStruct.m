function [varargout]=TestDeviationPopStruct(varargin)
% include_Optimisation
global TestDeviationPopStruct_Handle
try
nOut=nargout(TestDeviationPopStruct_Handle);
catch
include_Optimisation
nOut=nargout(TestDeviationPopStruct_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TestDeviationPopStruct_Handle(varargin{:});
end
