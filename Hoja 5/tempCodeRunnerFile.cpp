double t0 = 0; double tf = 10; int d = 100; string wname = "eval.dat";
    double t [d+1] {}; double y [d+1] {};
    eval(h, t0, tf, d, t, y);
    d_w_file_2cols(wname, t, y,d+1);
    //GnuplotPipe gp; gp.sendLine("plot '" + wname + "'");