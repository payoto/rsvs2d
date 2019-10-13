function [nameStr]=CallString2Name(callStr, mode)

	if nargin == 1
		mode = 'default';
	end

	switch mode
		case 'legacy'
			nameStr = regexprep(regexprep(regexprep(...
		        regexprep(callStr,'(\(|\)|,)','_'),'\.','_'),'''',''),...
		        ' ','');
		otherwise
			nameStr = matlab.lang.makeValidName(callStr);
	end

end