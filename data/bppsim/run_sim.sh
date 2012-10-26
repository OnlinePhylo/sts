for i in $(echo */); do
  cd $i
    echo $i
    bppseqgen param=bppseqgen.bpp
    # For phyml
    # seqmagick convert --pattern-replace '@' Z <in> <out>
    bppml param=bppml.bpp

    # phyml *.phyx 0 i 1 0 JC69 1.0 0.0 1 1.0 JC69.ML.dnd n n
  cd ..
done

grep -H likelihood */*.params.txt
