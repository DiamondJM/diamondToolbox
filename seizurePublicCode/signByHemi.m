function [arg1,arg2,arg3] = signByHemi(isLeftResection, arg1, arg2, arg3)

if isLeftResection
    arg1 = arg1 * -1; if nargin == 2; return; end 
    arg2 = arg2 * -1; if nargin == 3; return; end 
    arg3 = arg3 * -1; 
end 

end