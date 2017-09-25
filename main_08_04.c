//interface
//parse args, alloc mem, write/read to/from file
//must print overall solving time

/*
<path to bin file>
Usage: evc [input_file_name] [output_file_name] [options]
Where options include:
  -d    print debug messages [default OFF]
  -e    print errors [default OFF]
  -p    print matrix [default OFF]
  -t    print execution time [default OFF]
  -prec=<num>       precision [default - 1e-14]  
  -eps=<num>        'epsilon' [default - 1e-10]
  -max_iter=<num>   limit number of iterations 
                    [default - 0, i.e. not limit]
  -h, -?     print this and exit
*/