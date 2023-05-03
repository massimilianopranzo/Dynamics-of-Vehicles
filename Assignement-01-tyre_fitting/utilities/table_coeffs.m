function table_coeffs(tyre_coeffs)
    n = numel(fieldnames(tyre_coeffs));

    fid = fopen("table_coeffs.tex", 'w');
    cols = 6;
    rows = floor(n / (cols / 2));

    for i = 1: cols/2
        fprintf(fid, "Parameter & Value");
        if i ~= 3
            fprintf(fid, " & ");
        end
        if i == 3
            fprintf(fid, " \\\\\n");
        end
    end
    name = fieldnames(tyre_coeffs);
    
    k = 1;
    for j = 1:rows
        if name{k}(1) ~= 'L'
        fprintf(fid, "%s & ", name{k});
        fprintf(fid, "%.3f &", tyre_coeffs.(name{k}));
        k = k + 1;
        fprintf(fid, "%s & ", name{k});
        fprintf(fid, "%.3f &", tyre_coeffs.(name{k}));
        k = k + 1;
        fprintf(fid, "%s & ", name{k});
        fprintf(fid, "%.3f", tyre_coeffs.(name{k}));
        if mod(k, 3) == 0
            fprintf(fid, " \\\\\n");
        end
        k = k + 1;
        end
    end
    fclose(fid);
end