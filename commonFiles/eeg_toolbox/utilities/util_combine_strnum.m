function [combName] = util_combine_strnum(strArr, numArr)
numArr = arrayfun(@num2str, numArr, 'UniformOutput',false);
numArr(strcmpi(numArr, 'NaN')) = '';
combName = strcat(strArr, numArr);
return