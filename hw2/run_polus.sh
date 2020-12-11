mpisubmit.pl -p 10 -w 00:05 --stdout 10_128.txt bin/polus -- 1 1 1 128 128 128
mpisubmit.pl -p 20 -w 00:05 --stdout 20_128.txt bin/polus -- 1 1 1 128 128 128
mpisubmit.pl -p 40 -w 00:05 --stdout 40_128.txt bin/polus -- 1 1 1 128 128 128

mpisubmit.pl -p 10 -w 00:05 --stdout 10_256.txt bin/polus -- 1 1 1 256 256 256
mpisubmit.pl -p 20 -w 00:05 --stdout 20_256.txt bin/polus -- 1 1 1 256 256 256
mpisubmit.pl -p 40 -w 00:05 --stdout 40_256.txt bin/polus -- 1 1 1 256 256 256

mpisubmit.pl -p 10 -w 00:05 --stdout 10_512.txt bin/polus -- 1 1 1 512 512 512
mpisubmit.pl -p 20 -w 00:05 --stdout 20_512.txt bin/polus -- 1 1 1 512 512 512
mpisubmit.pl -p 40 -w 00:05 --stdout 40_512.txt bin/polus -- 1 1 1 512 512 512




mpisubmit.pl -p 10 -w 00:05 --stdout 10_128pi.txt bin/polus -- 3.14 3.14 3.14 128 128 128
mpisubmit.pl -p 20 -w 00:05 --stdout 20_128pi.txt bin/polus -- 3.14 3.14 3.14 128 128 128
mpisubmit.pl -p 40 -w 00:05 --stdout 40_128pi.txt bin/polus -- 3.14 3.14 3.14 128 128 128

mpisubmit.pl -p 10 -w 00:05 --stdout 10_256pi.txt bin/polus -- 3.14 3.14 3.14 256 256 256
mpisubmit.pl -p 20 -w 00:05 --stdout 20_256pi.txt bin/polus -- 3.14 3.14 3.14 256 256 256
mpisubmit.pl -p 40 -w 00:05 --stdout 40_256pi.txt bin/polus -- 3.14 3.14 3.14 256 256 256

mpisubmit.pl -p 10 -w 00:05 --stdout 10_512pi.txt bin/polus -- 3.14 3.14 3.14 512 512 512
mpisubmit.pl -p 20 -w 00:05 --stdout 20_512pi.txt bin/polus -- 3.14 3.14 3.14 512 512 512
mpisubmit.pl -p 40 -w 00:05 --stdout 40_512pi.txt bin/polus -- 3.14 3.14 3.14 512 512 512




