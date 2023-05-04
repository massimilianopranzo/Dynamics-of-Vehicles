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
            fprintf(fid, " \\\\ \\midrule \n");
        end
    end
    name = fieldnames(tyre_coeffs);
    
    k = 1;
    for j = 1:rows
        ncol = 1;
        while ncol <= cols/2
            if k > n
                break;
            else
                if all([(name{k} ~= "Fz01"), (name{k}(1) ~= 'L'), (name{k}(1) ~= 'Bz10')])
                    fprintf(fid, "%s & ", name{k});
                    if ncol == cols / 2
                        fprintf(fid, "%.4f  \t \\\\\n", tyre_coeffs.(name{k}));
                    else
                        fprintf(fid, "%.4f & ", tyre_coeffs.(name{k}));
                    end
                    ncol = ncol + 1;
                    
                    k = k + 1;
                else
                    k = k + 1;
                end
            end
        end
    end
    fprintf(fid, "\\bottomrule");
    fclose(fid);
end