function output = madTK(input)
% function output = madTK(input, xnan)
% xnan, 1=exclude nan, 2=include nan

% if nargin == 1;
%     xnan = 0;
% end
% 
% if size(input,1) > 1 && size(input,2)
%     error;
% elseif size(input,1) > 1 
%     input = input';
% end
% 
% switch xnan
%     case 0
%         output = 1.4826*(median(abs(input - median(input))));
%     case 1
%         output = 1.4826*(nanmedian(abs(input - nanmedian(input))));
% end
warndlg('this function is broken');
output = 2.5*1.4826*(nanmedian(abs(input - nanmedian(input))));