function [broken_ind]=findgaps(input_ind)
% finds gaps (greater than 0) in between indeces in a vector

gaps =find(diff(input_ind) > 1);
broken_ind = zeros(length(gaps),50);
if ~isempty(gaps)
    for gapind = 1:length(gaps)+1;
        if gapind == 1;
            temp = input_ind(1:gaps(gapind));
        elseif gapind == length(gaps)+1
            temp = input_ind(gaps(gapind-1)+1:end);
        else
            temp = input_ind(gaps(gapind-1)+1:gaps(gapind));
        end
        broken_ind(gapind,1:length(temp)) = temp;
    end
end

end