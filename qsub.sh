 qsub -cwd -l vf=2g -P st_pg -q st.q *.sh
 qstat
 qdel
