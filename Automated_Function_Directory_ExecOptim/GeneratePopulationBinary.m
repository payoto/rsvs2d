function [varargout]=GeneratePopulationBinary(varargin)
% OptimisationOutput
global GeneratePopulationBinary_Handle
nOut=nargout(GeneratePopulationBinary_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GeneratePopulationBinary_Handle(varargin{:});
end
