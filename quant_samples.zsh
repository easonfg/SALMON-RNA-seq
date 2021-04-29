
run_salmon(){
  echo "Processing sample $1"
  salmon quant --gcBias -i /Users/hhuang/Desktop/ketone_exp/salmon/mus_musculus_index -l A \
           -1 data/$1R1_001.fastq.gz \
           -2 data/$1R2_001.fastq.gz \
           -p 16 \
           --validateMappings \
            --useVBOpt \
            --seqBias \
           -o quants/$1_quant
}


#for counter in {0..4}
#do
#  i=$((counter * 2))
#  #echo $i
#  run_salmon KD_P${i}_1
#  run_salmon KD_P${i}_2
#  run_salmon S_P${i}_1
#  run_salmon S_P${i}_2
#done


for counter in A1_339_S304_L004_ A2_340_S305_L004_ A3_341_S306_L004_ A4_338_S307_L004_ C1_322_S300_L004_ C2_316_S301_L004_ C3_356_S302_L004_ C4_355_S303_L004_
do
  run_salmon ${counter}
done
