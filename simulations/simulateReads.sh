


genename=TRA


# imrep_revision2/$genename/trans.fa is a fasta file with the immune transcripts


for k in 1 2 4 8 16 32 64 128; do
    for l in 50 75 100; do
        ins=$(( 3 * $l ))
        ./simNGS/bin/simLibrary -r $l -i $ins -x $k imrep_revision2/$genename/trans.fa > imrep_revision2/$genename/ref.frag
        ./simNGS/bin/simNGS runfiles/$l.runfile imrep_revision2/$genename/ref.frag -p paired > imrep_revision2/$genename/sim_rl_${l}_cov_${k}.fq
    done
done





