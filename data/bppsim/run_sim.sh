for i in $(echo */); do
  cd $i
    echo $i
    bppseqgen param=bppseqgen.bpp
    bppml param=bppml.bpp
  cd ..
done

grep -H likelihood */*.params.txt
