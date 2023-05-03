function table_P0(action, name_coeffs, P0, lb, ub)
    cols = max(cellfun(@numel, P0));
    rowsetts = 4; % Rows for each minimization: P0, lb, ub
    n_minim = length(P0);
    name_arrays = {"Name","Guess", "Lower Bound", "Upper Bound"};
    fid = fopen(strcat("P0_",action,".tex"), 'w');
    fprintf(fid, "\\begin{tabular}{%s|" + "%s" + "} \n", 'c', repmat('c', cols, 1));
    fprintf(fid, "\\toprule \n");
    idx = 1;
    for i = 1:n_minim
        for k = 1:cols
            if k == 1
                fprintf(fid, "%s & ", name_arrays{1});
            end
            if k > length(name_coeffs{idx})
                if k == cols
                    fprintf(fid, "\\\\ \n");
                else
                    fprintf(fid, "& ");
                end
            else
                if k ~= cols
                    fprintf(fid, "%s & ", name_coeffs{idx}(k));
                else
                    fprintf(fid, "%s \\\\ \n", name_coeffs{idx}(k));
                end
            end
        end
        for k = 1:cols
            if k == 1
                fprintf(fid, "%s & ", name_arrays{2});
            end
            if k > length(P0{idx})
                fprintf(fid, "& ");
                if k == cols
                    fprintf(fid, "\\\\ \n");
                end
            else
                if k ~= cols
                    fprintf(fid, "%.3f & ", P0{idx}(k));
                else
                    fprintf(fid, "%.3f \\\\ \n", P0{idx}(k));
                end
            end
        end
        for k = 1:cols
            if k == 1
                fprintf(fid, "%s & ", name_arrays{3});
            end
            if isempty(lb{idx}) || k > length(lb{idx})
                fprintf(fid, "& ");
                if k == cols
                    fprintf(fid, "\\\\ \n");
                end
            else
                if k ~= cols
                    if isinf(lb{idx}(k))
                        fprintf(fid, "$-\\infty$ & ");
                    else
                        fprintf(fid, "%.3f & ", lb{idx}(k));
                    end
                else
                    if isinf(lb{idx}(k))
                        fprintf(fid, "$-\\infty$ \\\\ \n");
                    else
                        fprintf(fid, "%.3f \\\\ \n", lb{idx}(k));
                    end
                end
            end
        end
        for k = 1:cols
            if k == 1
                fprintf(fid, "%s & ", name_arrays{4});
            end
            if isempty(ub{idx}) || k > length(ub{idx})
                fprintf(fid, "& ");
                if k == cols
                    fprintf(fid, "\\\\ \n");
                end
            else
                if k ~= cols
                    if isinf(ub{idx}(k))
                        fprintf(fid, "$\\infty$ & ");
                    else
                        fprintf(fid, "%.3f & ", ub{idx}(k));
                    end
                else
                    if isinf(ub{idx}(k))
                        fprintf(fid, "$\\infty$ \\\\ \n");
                    else
                        fprintf(fid, "%.3f \\\\ \n", ub{idx}(k));
                    end
                end
            end
        end
        fprintf(fid, "\\midrule \n");
        idx = idx + 1;
    end    
    fprintf(fid, "\\bottomrule \n");
    fprintf(fid, "\\end{tabular} \n");
    fclose(fid);
end
