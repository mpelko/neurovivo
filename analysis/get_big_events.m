arg_list = argv ();

load(arg_list{1});
cd(analysis_folder);
be = bigevents(min_amplitude, timelimit, 1000*voltage, times);
amplitudes = be.ampli;
start_bins = be.binstart;
result = [amplitudes; start_bins];
save("-7", arg_list{1},"result");
%be.ampli;