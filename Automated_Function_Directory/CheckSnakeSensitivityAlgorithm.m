function [varargout]=CheckSnakeSensitivityAlgorithm(varargin)
% include_Optimisation
global CheckSnakeSensitivityAlgorithm_Handle
try
nOut=nargout(CheckSnakeSensitivityAlgorithm_Handle);
catch
include_Optimisation
nOut=nargout(CheckSnakeSensitivityAlgorithm_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CheckSnakeSensitivityAlgorithm_Handle(varargin{:});
end
