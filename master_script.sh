#!/bin/bash

# this script is used to write individual running scripts and steps.

YJI_SPACE=/run/media/qu/7de37f23-d1d5-4201-b381-28d1d2810570/yji/

cd $YJI_SPACE/bgi/01
ls -d */ | sed 's/\///' > all_folders_01.list

mkdir raw_fq
mkdir trimmed_reads
mkdir B10K_fa
mkdir fa
mkdir done_trimmed-reads

while read folder; do
    echo '#!/bin/bash' > 0_trimmomatic_${folder}.sh
    echo "cd $folder" >> 0_trimmomatic_${folder}.sh
    echo 'for file in `ls *.1.fq.gz`; do
        base=${file%%.1.fq.gz}
        gunzip -c $file  > ../raw_fq/${base}.1.cp.fq
        gunzip -c ${base}.2.fq.gz > ../raw_fq/${base}.2.cp.fq
        cd .. ' >> 0_trimmomatic_${folder}.sh
    echo '
        java -jar /run/media/qu/7de37f23-d1d5-4201-b381-28d1d2810570/share/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
            PE -threads 4 -phred33 \
            raw_fq/${base}.1.cp.fq \
            raw_fq/${base}.2.cp.fq \
            trimmed_reads/${base}_p1.fq.gz trimmed_reads/${base}_u1.fq.gz \
            trimmed_reads/${base}_p2.fq.gz trimmed_reads/${base}_u2.fq.gz \
            LEADING:20 TRAILING:20 MINLEN:10' >>  0_trimmomatic_${folder}.sh
    echo '
        echo $base trimmomatic done' >> 0_trimmomatic_${folder}.sh
    echo 'done' >> 0_trimmomatic_${folder}.sh
    
    echo "cd ../trimmed_reads" >> 0_trimmomatic_${folder}.sh
    echo 'speciesNum=0' >> 0_trimmomatic_${folder}.sh
    echo '
for file in `ls *.fq.gz`; do
              base=${file%%.fq.gz}
              id=${base%%_*}
              species=`grep $id ../id-species.txt | cut -f2`
              gunzip $file
              paste - - - - < ${base}.fq | cut -f1,2 | \
              sed ' "'" "s/^@/>/" "'" '| tr "\t" "\n" | sed ' "'" "s/ /-/; s/#/-/" "'" ' >  ../B10K_fa/${species}_${base}.fa' >> 0_trimmomatic_${folder}.sh
    echo '    baseNum=`grep -c ">" ../B10K_fa/${species}_${base}.fa` ' >> 0_trimmomatic_${folder}.sh
    echo '    ((speciesNum=speciesNum+baseNum))' >> 0_trimmomatic_${folder}.sh
    echo 'done' >> 0_trimmomatic_${folder}.sh

    echo 'cd ..' >> 0_trimmomatic_${folder}.sh
    echo 'cat B10K_fa/${species}_*.fa > fa/${species}.fa' >> 0_trimmomatic_${folder}.sh

    echo 'totalNum=`grep -c ">" fa/${species}.fa`' >> 0_trimmomatic_${folder}.sh
    echo 'if [ "$totalNum" -eq "$speciesNum" ]; then ' >> 0_trimmomatic_${folder}.sh
    echo '    echo "sequence number matched for $species"' >> 0_trimmomatic_${folder}.sh
    echo '    rm raw_fq/*.fq' >> 0_trimmomatic_${folder}.sh
    echo '    rm B10K_fa/*.fa' >> 0_trimmomatic_${folder}.sh
    echo '    mv trimmed_reads/*.fq done_trimmed-reads/' >> 0_trimmomatic_${folder}.sh
    echo 'else' >> 0_trimmomatic_${folder}.sh
    echo '    echo "sequence number doesnt match for $species"' >> 0_trimmomatic_${folder}.sh
    echo '    exit 1' >> 0_trimmomatic_${folder}.sh
    echo 'fi' >> 0_trimmomatic_${folder}.sh
done < all_folders_01.list 

cp 0_trimmomatic_*.sh ../trimmomatic_scripts
cd trimmomatic_scripts
i=1
for job in 0_trimmomatic_*.sh; do
#    if [ $i -gt 10 -a $i -lt 21 ]; then
    if [ $i -gt 20 -a $i -lt 31 ]; then # changed 4.23
        echo $job
         cp $job ../
        ((i=i+1))
    else
        ((i=i+1))
    fi
done

cd ..
for job in 0_trimmomatic_*.sh; do
    bash $job > ${job}.out 2>&1 &
done

cd done_trimmed-reads
for fq in *.fq; do
     gzip $fq 
    echo "$fq gzip trimmed fq done"
done
cd ..

mkdir trf
cd trf
for fasta in `ls ../fa/*.fa`; do
    species0=${fasta%%.fa}
    species=${species0##*/}
    trf_base=${species}.fa.2.5.7.80.10.20.7 
    echo '#!/bin/bash' > ${species}_trf.sh
    echo "trf $fasta 2 5 7 80 10 20 7 -f -h > /dev/null 2>&1" >> ${species}_trf.sh
    #to be run
    echo 'dat_nseq=`grep -c Sequence: ' "${species}.fa.2.5.7.80.10.20.7.dat" '`' >> ${species}_trf.sh
    echo 'nseq=`grep -c ">"' " $fasta " '`' >> ${species}_trf.sh 
    echo 'if [ $nseq -eq $dat_nseq ]; then' >> ${species}_trf.sh
    echo '    echo number of sequence match for '"$species" >> ${species}_trf.sh
    echo "    trf-sort.pl ${trf_base}.dat > ${trf_base}.sorted.dat" >> ${species}_trf.sh
    echo "    trf-redundancy_v2.pl -i ${trf_base}.sorted.dat -f" >> ${species}_trf.sh
    echo '    awk '"'"'$5>=5'"'"' \' >> ${species}_trf.sh
    echo "    ${trf_base}.sorted.mdat > ${trf_base}.sorted.als5.mdat" >> ${species}_trf.sh
    
    echo "    column-transfer.pl -i ${trf_base}.sorted.als5.mdat,15,17 \
                  -r $YJI_SPACE/bgi/msat1-7_L110_collapsed.txt,3,2 \
                  > ${trf_base}.sorted.als5.c17.mdat " >> ${species}_trf.sh
    echo "    column-transfer.pl -i ${trf_base}.sorted.als5.c17.mdat,17,18 \
                  -r $YJI_SPACE/bgi/msat1-7_L110_collapsed.txt,3,2 \
                  > ${trf_base}.sorted.als5.c18.mdat " >> ${species}_trf.sh

    echo "    rm ${trf_base}.sorted.dat" >> ${species}_trf.sh
    echo "    rm ${trf_base}.sorted.mdat" >> ${species}_trf.sh
    echo "    rm ${trf_base}.sorted.als5.mdat" >> ${species}_trf.sh
    echo "    rm ${trf_base}.sorted.als5.c17.mdat" >> ${species}_trf.sh
    echo "    nohup gzip ${trf_base}.dat &" >> ${species}_trf.sh
    echo "fi" >> ${species}_trf.sh
done

## 1. prepare 15-bp flanks
mkdir 1_prepare-flank-15bp
cd 1_prepare-flank-15bp
while read motif; do
    mkdir $motif
done < ../../all_motif.list

for file in `ls ../trf/*.sorted.als5.c18.mdat`; do
    file_name=${file##*/}
    species=${file_name%%.fa.2.5.7.80.10.20.7.sorted.als5.c18.mdat}
    echo '#!/bin/bash' > ${species}_1-prepare-flank-15bp.sh
    while read motif; do
        echo "cd $motif" >> ${species}_1-prepare-flank-15bp.sh
        echo "awk -v motif="'"'$motif'"' "'" '$18==motif' "'" "../../trf/${species}.fa.2.5.7.80.10.20.7.sorted.als5.c18.mdat \
         > ${species}_${motif}-only_re-fmt.sorted.als5.c18.mdat" >> ${species}_1-prepare-flank-15bp.sh
        echo "cut -f1 ${species}_${motif}-only_re-fmt.sorted.als5.c18.mdat | sort | uniq > ${species}_${motif}-only.list" >> ${species}_1-prepare-flank-15bp.sh 
        echo "faSomeRecords ../../fa/${species}.fa ${species}_${motif}-only.list ${species}_${motif}-only.fasta" >> ${species}_1-prepare-flank-15bp.sh
        echo "bioawk -c fastx '" '{print $name, length($seq) }' "'" "< ${species}_${motif}-only.fasta > ${species}_${motif}-only_length.txt" >> ${species}_1-prepare-flank-15bp.sh
        echo "mdat_flank_filter.pl -i ${species}_${motif}-only_re-fmt.sorted.als5.c18.mdat \
                -l 15 -t ${species}_${motif}-only_length.txt \
                -o ${species}_${motif}-only " >> ${species}_1-prepare-flank-15bp.sh
        echo "mdat_flank2bed.pl -i ${species}_${motif}-only_w-both-flanks.mdat -l 15 \
                -t ${species}_${motif}-only_length.txt \
                > ${species}_${motif}-only_both-flanks.bed" >> ${species}_1-prepare-flank-15bp.sh
        echo "clip_bed.pl -i ${species}_${motif}-only_both-flanks.bed \
                 -s ${species}_${motif}-only.fasta \
                 > ${species}_${motif}-only_both-flanks.fasta" >> ${species}_1-prepare-flank-15bp.sh
        echo "cat_flanks.pl ${species}_${motif}-only_both-flanks.fasta \
                 > ${species}_${motif}-only_cat-flanks.fasta" >> ${species}_1-prepare-flank-15bp.sh
        echo "rm ${species}_${motif}-only.list " >> ${species}_1-prepare-flank-15bp.sh
        echo "rm ${species}_${motif}-only_both-flanks.bed " >> ${species}_1-prepare-flank-15bp.sh
        echo "rm ${species}_${motif}-only_both-flanks.fasta " >> ${species}_1-prepare-flank-15bp.sh
        echo "cd .." >> ${species}_1-prepare-flank-15bp.sh
        echo "echo $species $motif done" >>  ${species}_1-prepare-flank-15bp.sh
    done < ../../all_motif.list
done 

mkdir 2_self-blast
cd 2_self-blast
while read motif; do
    mkdir $motif
done < ../../all_motif.list
          
for file in `ls ../trf/*.sorted.als5.c18.mdat`; do
    file_name=${file##*/}
    species=${file_name%%.fa.2.5.7.80.10.20.7.sorted.als5.c18.mdat}
    echo '#!/bin/bash' > ${species}_2-self-blast.sh
    while read motif; do
        echo "cd $motif" >> ${species}_2-self-blast.sh
        echo "ln -s ../../1_prepare-flank-15bp/${motif}/${species}_${motif}-only_cat-flanks.fasta ." >> ${species}_2-self-blast.sh
        echo "filterN.pl -s ${species}_${motif}-only_cat-flanks.fasta -n 2" >> ${species}_2-self-blast.sh 
        echo "makeblastdb -in ${species}_${motif}-only_cat-flanks_filtered-N.fasta -dbtype nucl" >> ${species}_2-self-blast.sh
        echo "blastn -task blastn-short \
        -query ${species}_${motif}-only_cat-flanks_filtered-N.fasta \
        -db ${species}_${motif}-only_cat-flanks_filtered-N.fasta \
        -outfmt 6 \
        -qcov_hsp_perc 100 \
        -out ${species}_${motif}-only_cat-flanks_all-v-all.fmt6 " >> ${species}_2-self-blast.sh
        echo "#silix -i 0.93 -r 0.99 -l 29 -m 0.99 \
              #    ${species}_${motif}-only_cat-flanks_filtered-N.fasta \
              #    ${species}_${motif}-only_cat-flanks_all-v-all.fmt6 \
              #    > ${species}_${motif}-only_cat-flanks_silix.fnodes" >> ${species}_2-self-blast.sh
        echo "#mkdir ${species}_${motif}-only_silix_clusters" >> ${species}_2-self-blast.sh
        echo "#cd ${species}_${motif}-only_silix_clusters" >> ${species}_2-self-blast.sh
        echo "#silix-split ../${species}_${motif}-only_cat-flanks_filtered-N.fasta \
              #    ../${species}_${motif}-only_cat-flanks_silix.fnodes" >> ${species}_2-self-blast.sh
        echo "#for i in *.fasta; do " >> ${species}_2-self-blast.sh
        echo '#  num=`grep -c ">" $i`' >> ${species}_2-self-blast.sh
        echo '#  if [ "$num" -le 100 ]; then ' >> ${species}_2-self-blast.sh
        echo '#      head -n2 $i >> ' "../${species}_${motif}-only_cat-flanks_reduced-100.fasta" >> ${species}_2-self-blast.sh
        echo '#  fi' >> ${species}_2-self-blast.sh
        echo '#done' >> ${species}_2-self-blast.sh
        echo "#cd .." >> ${species}_2-self-blast.sh
        echo "#tar -cf ${species}_${motif}-only_silix_clusters.tar ${species}_${motif}-only_silix_clusters" >> ${species}_2-self-blast.sh
        echo '#if [ $? -eq 0 ]; ' "then rm -r ${species}_${motif}-only_silix_clusters; fi" >> ${species}_2-self-blast.sh 
        echo "#rm ${species}_${motif}-only_cat-flanks_filtered-N.fasta.*" >> ${species}_2-self-blast.sh
        echo "#rm ${species}_${motif}-only_cat-flanks.fasta" >> ${species}_2-self-blast.sh
        echo "cd .." >> ${species}_2-self-blast.sh
        echo "echo $species $motif self blast done" >>  ${species}_2-self-blast.sh
    done < ../../all_motif.list
done

mkdir 3_blast-galGal
cd 3_blast-galGal
while read motif; do
    mkdir $motif
done < ../../all_motif.list

while read motif; do
    echo '#!/bin/bash' > ${motif}_3-blast-galGal.sh
    echo "cd $motif" >> ${motif}_3-blast-galGal.sh
    for file in `ls ../2_self-blast/${motif}/*_reduced-100.fasta`; do
        filename=${file##*/}
        species=${filename%_*-only_cat-flanks_reduced-100.fasta}
        echo "ln -s ../../2_self-blast/${motif}/$filename ." >> ${motif}_3-blast-galGal.sh
        echo "makeblastdb -dbtype nucl -in $filename " >> ${motif}_3-blast-galGal.sh
        echo "blastn -task blastn \
                  -query $YJI_SPACE/bgi/galGal-${motif}_loci-20.fasta \
                  -db $filename \
                  -word_size 4 -reward 1 -penalty -1 \
                  -qcov_hsp_perc 100 \
                  -outfmt 6 \
                  -out galGal_${species}_${motif}_blastn.fmt6 " >> ${motif}_3-blast-galGal.sh
        echo "cut -f1  galGal_${species}_${motif}_blastn.fmt6| uniq -u > galGal_${species}_${motif}_blastn_unique-query.list" >> ${motif}_3-blast-galGal.sh 
        echo "grep -f  galGal_${species}_${motif}_blastn_unique-query.list \
                 galGal_${species}_${motif}_blastn.fmt6 \
                > galGal_${species}_${motif}_blastn_unique-query.fmt6" >> ${motif}_3-blast-galGal.sh
        echo 'fmt6_num=`wc -l '  "galGal_${species}_${motif}_blastn.fmt6 | cut -f1 -d" '" "`' >>  ${motif}_3-blast-galGal.sh
        echo 'uniqQuery_num=`wc -l ' " galGal_${species}_${motif}_blastn_unique-query.fmt6 | cut -f1 -d" '" "`'  >>  ${motif}_3-blast-galGal.sh
        echo 'if [[ $fmt6_num -ne $uniqQuery_num ]]; then' >> ${motif}_3-blast-galGal.sh 
        echo "    echo $species $motif has redundant queries of galGal " >>  ${motif}_3-blast-galGal.sh
        echo "else" >>  ${motif}_3-blast-galGal.sh
        echo "    rm galGal_${species}_${motif}_blastn.fmt6 " >> ${motif}_3-blast-galGal.sh
        echo "fi" >> ${motif}_3-blast-galGal.sh
        echo "rm $filename " >> ${motif}_3-blast-galGal.sh
        echo "rm ${filename}.* " >> ${motif}_3-blast-galGal.sh 
        echo "rm galGal_${species}_${motif}_blastn_unique-query.list " >> ${motif}_3-blast-galGal.sh
        echo "echo ${species} ${motif} blast galGal done" >> ${motif}_3-blast-galGal.sh
     done 
     echo "cd .." >> ${motif}_3-blast-galGal.sh
done < ../../all_motif.list

mkdir 4_genotype
cd 4_genotype
while read motif; do
    mkdir $motif
done < ../../all_motif.list


#loop through each motif, each galGal-locus and each species.
while read motif; do
    cd $motif
    while read galGal_locus; do
        locusID=`grep $galGal_locus ../../../all_loci_info.txt | cut -f2 `
        for file in `ls ../../3_blast-galGal/${motif}/*_unique-query.fmt6`; do
            filename=${file##*/}
            species=`cut -f2,3 -d"_" <<< $filename`
            if [ "$(grep -c $galGal_locus $file)" -ne 0 ]; then
                species_locus=($(grep $galGal_locus $file | cut -f2))
                cluster_id=`awk -v pattern="$species_locus" '$2~ pattern' ../../2_self-blast/${motif}/${species}_${motif}-only_cat-flanks_silix.fnodes | cut -f1`
                if [ ! -d "$species_${motif}-only_silix_clusters" ]; then
                    tar -xf ../../2_self-blast/${motif}/${species}_${motif}-only_silix_clusters.tar
                fi
                cd ${species}_${motif}-only_silix_clusters
                grep ">" ${species}_${motif}-only_cat-flanks_filtered-N_${cluster_id}.fasta |\
                    sed 's/>//' | cut -f1,2,5 -d"_" | awk -F"_" '{print $1 "\t" $2-1 "\t" $3}' \
                    > ${species}_${motif}-only_cat-flanks_filtered-N_${cluster_id}.bed
                cut -f1 ${species}_${motif}-only_cat-flanks_filtered-N_${cluster_id}.bed \
                    > ${species}_${motif}-only_${cluster_id}_prefa.list
                faSomeRecords ../../../1_prepare-flank-15bp/${motif}/${species}_${motif}-only.fasta \
                    ${species}_${motif}-only_${cluster_id}_prefa.list \
                    ${species}_${motif}-only_${cluster_id}.fasta
                clip_bed.pl \
                    -i ${species}_${motif}-only_cat-flanks_filtered-N_${cluster_id}.bed \
                    -s ${species}_${motif}-only_${cluster_id}.fasta \
                    > ../${species}_${locusID}.fasta
                cd ..
                sed -i "s/>/>${species}_${locusID}_/" ${species}_$locusID.fasta
                echo $species $locusID of $motif extract sequence done
            fi
        done
    done < ../../../galGal-${motif}_loci-20.list
    rm -r *_${motif}-only_silix_clusters
    cd ..
#done < ../../all_motif_2_self-blast-done.list
done < ../../AA-AT_motif.list

cat /dev/null > all-motifs_clipped_length.txt
while read motif; do
    cd $motif
    while read galGal_locus; do
        locusID=`grep $galGal_locus ../../../all_loci_info.txt | cut -f2 `
        cat *_${locusID}.fasta $YJI_SPACE/bgi/all_species_data/${locusID}_${motif}_all-sp.fasta > ${locusID}_${motif}_combined.fasta
        mafft --adjustdirection ${locusID}_${motif}_combined.fasta > ${locusID}_${motif}_combined_mafft.fasta
        python3 $YJI_SPACE/python/clip_msa.py \
        -i ${locusID}_${motif}_combined_mafft.fasta \
        -n 15 \
        > ${locusID}_${motif}_combined_mafft_clipped.fasta
        sed -i 's/:/./g; s/\//./; s/_/-/g' ${locusID}_${motif}_combined_mafft_clipped.fasta 
        perl $YJI_SPACE/perl/getLength.pl ${locusID}_${motif}_combined_mafft_clipped.fasta \
        >> ../all-motifs_clipped_length.txt
        echo $motif $locus done
        rm *_${locusID}.fasta
        rm ${locusID}_${motif}_combined_mafft.fasta
        rm ${locusID}_${motif}_combined_mafft_clipped.fasta.index
    done < ../../../galGal-${motif}_loci-20.list
    cd ..
done < ../../AA-AT_motif.list


perl $YJI_SPACE/perl/reformat_seqLength2STR-FM.pl \
    all-motifs_clipped_length.txt ../../all_loci_info.txt 5 \
    > all-motifs_STR-FM-in_d5.txt

python $YJI_SPACE/python/STR-FM/GenotypeTRcorrection_yj.py \
    all-motifs_STR-FM-in_d5.txt \
    $YJI_SPACE/python/STR-FM/PCRinclude.allrate.bymajorallele_modified.txt \
    all-motifs_STR-FM-out_d5.txt \
    0.5
sed -i 's/-L/_L/' all-motifs_STR-FM-out_d5.txt 

python ~/lib/STR-FM/GenotypeTRcorrection_yj.py \
    all-motifs_STR-FM-in_d4.txt \
    ~/lib/STR-FM/PCRinclude.allrate.bymajorallele_modified.txt \
    all-motifs_STR-FM-out_d4.txt \
    0.5


