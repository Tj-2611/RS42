function [model, filename] = find_model(year, hy, flag)
    if strcmp(flag, '42rswithmu')
        if year == 2020
            filename = [num2str(year), '_half', num2str(hy), '_', flag, '.mat'];
            if year == 2020 && hy == 2
                model = '42rswithmu';
            else
                model = '42withmu'; 
            end
        else
            filename = [num2str(year), '_', flag, '.mat'];
            if ismember(year, [2015, 2017, 2019])
                model = '42withmu';
            else
                model = '42rswithmu';
            end
        end
    else
        model = flag;
        if year == 2020
            filename = [num2str(year), '_half', num2str(hy), '_', flag, '.mat'];
        else
            filename = [num2str(year), '_', flag, '.mat'];
        end
    end

end

