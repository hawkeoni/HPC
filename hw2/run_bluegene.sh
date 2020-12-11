mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 128 128 128 --stdout no_omp_out_proc_64_grid_128.txt
mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 128 128 128 --stdout no_omp_out_proc_128_grid_128.txt
mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 128 128 128 --stdout no_omp_out_proc_256_grid_128.txt


mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 256 256 256 --stdout no_omp_out_proc_64_grid_256.txt
mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 256 256 256 --stdout no_omp_out_proc_128_grid_256.txt
mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 256 256 256 --stdout no_omp_out_proc_256_grid_256.txt

mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 512 512 512 --stdout no_omp_out_proc_64_grid_512.txt
mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 512 512 512 --stdout no_omp_out_proc_128_grid_512.txt
mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=1" bin/bluegene 1 1 1 512 512 512 --stdout no_omp_out_proc_256_grid_512.txt



mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 128 128 128 --stdout omp_out_proc_64_grid_128.txt
mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 128 128 128 --stdout omp_out_proc_128_grid_128.txt
mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 128 128 128 --stdout omp_out_proc_256_grid_128.txt


mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 256 256 256 --stdout omp_out_proc_64_grid_256.txt
mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 256 256 256 --stdout omp_out_proc_128_grid_256.txt
mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 256 256 256 --stdout omp_out_proc_256_grid_256.txt

mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 512 512 512 --stdout omp_out_proc_64_grid_512.txt
mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 512 512 512 --stdout omp_out_proc_128_grid_512.txt
mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=4" bin/bluegene 1 1 1 512 512 512 --stdout omp_out_proc_256_grid_512.txt

#mpisubmit.bg -n 1 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 128 128 128 --stdout no_omp_out_proc_1_grid_128.txt
#mpisubmit.bg -n 2 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 128 128 128 --stdout no_omp_out_proc_2_grid_128.txt
#mpisubmit.bg -n 4 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 128 128 128 --stdout no_omp_out_proc_4_grid_128.txt
#mpisubmit.bg -n 8 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 128 128 128 --stdout no_omp_out_proc_8_grid_128.txt
#
#mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 256 256 256 --stdout no_omp_out_proc_64_grid_256.txt
#mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 256 256 256 --stdout no_omp_out_proc_128_grid_256.txt
#mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 256 256 256 --stdout no_omp_out_proc_256_grid_256.txt
#
#mpisubmit.bg -n 64 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 512 512 512 --stdout no_omp_out_proc_64_grid_512.txt
#mpisubmit.bg -n 128 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 512 512 512 --stdout no_omp_out_proc_128_grid_512.txt
#mpisubmit.bg -n 256 -w 00:10:00 -e "OMP_NUM_THREADS=2" bin/bluegene 3.14 3.14 3.14 512 512 512 --stdout no_omp_out_proc_256_grid_512.txt
