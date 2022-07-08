args = argv();
%printf("args = %s\n", args{1});

function z = justrandomfunction(func_args)
	k = str2double(func_args);
    z = k^.5 + 1.5;
end

for i = 1:5
	out{i} = justrandomfunction(args{i});
end

disp(out)
%x = -10:0.1:10;
%plot (x, sin (x));
%pause(5)