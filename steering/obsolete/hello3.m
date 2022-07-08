args = argv();
printf("args = %s\n", args{1});

function z = justrandomfunction(func_args)
	k = str2double(func_args);
    z = k^.5 + 1.5
end

out = justrandomfunction(args{1})

%x = -10:0.1:10;
%plot (x, sin (x));
%pause(5)