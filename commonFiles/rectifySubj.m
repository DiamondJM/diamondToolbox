function [twoDigit,sixDigit] = rectifySubj(subj)


if isa(subj,'cell')
    sixDigit = cell(size(subj));
    twoDigit = zeros(size(subj)); 
    for ii = 1:numel(subj)
        [twoDigit(ii),sixDigit{ii}] = rectifySubj(subj{ii});
    end
    return
end 

if ischar(subj)
    if length(subj) == 2
        sixDigit = sprintf('NIH0%s',subj);
        twoDigit = str2double(subj);
    else
        sixDigit = subj; 
        twoDigit = str2double(subj(end-1:end)); 
    end
else
    twoDigit = subj; 
    sixDigit = sprintf('%3d',twoDigit); 
    sixDigit(ismember(sixDigit,' ')) = '0'; 

    sixDigit = sprintf('NIH%s',sixDigit);
end

end
