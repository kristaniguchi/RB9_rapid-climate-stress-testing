### Written by Donny Kim (dkim8398@sdsu.edu or donnykim27@gmail.com) in June, 2022
### Gegraphy Department, San Diego State University

function output = preprocess_table2csv(table, TS_name)
	% converting table to cell, and then appending TS_name and VariableNames(colnames) on top two rows of the cell.
	x = table2cell(table);
	y = table.Properties.VariableNames;
	y2 = [repmat({TS_name}, [1, size(y, 2)])];
	output = vertcat(y2, y, x);
end