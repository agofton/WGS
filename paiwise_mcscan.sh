#!/bin/bash

# All output will be in the current directory

# List genomes here. Each should have a .gff and .fna and be in the current dir.
g1=AFGY1_genomic
g2=afzelii_Pko_genomic
g3=anserina_BA2_genomic
g4=bissettiae_DN127genomic
g5=Btach_genomic
g6=Bturcica_genomic
g7=burgdorferi_BB31_genomic
g8=coriaceae_Co53_genomic
g9=duttonii_Ly_genomic
g10=finlandensis_SV1_genomic
g11=garinii_20047_genomic
g12=hermsi_HS1_genomic
g13=mayonii_MN14-1539_genomic
g14=miyamotoi_Izh-4_genomic
g15=parkeri_SLO_genomic
g16=recurrentis_A1_genomic
g17=RT1S_genomic
g18=RT5S_genomic
g19=spielmanii_A14S_genomic
g20=turicatae_91E135_genomic
g21=valaisiana_ABCY02_genomic

list=( "$g1" "$g2" "$g3" "$g4" "$g5" "$g6" "$g7" "$g8" "$g9" "$g10" "$g11" "$g12" "$g13" "$g14" "$g15" "$g16" "$g17" "$g18"
	"$g19" "$g20" "$g21" )

# For each genome file, get chromosome names and lengths
for x in ${list[@]}; do
	echo "Creating genome.length files" && echo ""
	python -m jcvi.formats.gff bed --type=region --key=Name ${x}.gff | cut -f1,3 | sed -e 's/\s\+/,/g' > ${x}.genome.lengths
done
# concatenate genome lenghts files and create headers
echo 'chromosome,size' > all.genome.lengths
cat *_genomic.genome.lengths >> all.genome.lengths

# Reformat .gff as .bed
for x in ${list[@]}; do
	echo "Reformatting .gff as .bed" && echo ""
	python -m jcvi.formats.gff bed --type=CDS --key=Name ${x}.gff -o ${x}.bed
done

# Extract CDS from .fna
for x in ${list[@]}; do
	echo "Extracting CDS from .bed" && echo ""
	bedtools getfasta -fi ${x}.fna -fo ${x}.cds -bed ${x}.bed -nameOnly
done

# Start pairwise synteny searches
for x in ${list[@]}; do
	for y in ${list[@]}; do
		echo "Performing ortholog search" && echo ""
		python -m jcvi.compara.catalog ortholog $x $y --no_strip_names --cscore=.99
		echo "Scanning for syntenic blocks" && echo ""
		python -m jcvi.compara.synteny screen --simple ${x}.${y}.anchors ${x}.${y}.anchors.new

		# Convert Genes to chromosomes & bp coordinates
		anch=${x}.${y}.last.filtered

		# Get Gene 1 coords
		cut -f1 ${anch} > ${anch}.g1

		while read line; do
			grep -m 1 ${line} ${x}.bed | cut -f1,2,3 >> ${anch}.g1.coords
		done<${anch}.g1

		# Get gene 2 coords
		cut -f2 ${anch} > ${anch}.g2

		while read line; do
			grep -m 1 ${line} ${y}.bed | cut -f1,2,3 >> ${anch}.g2.coords
		done<${anch}.g2

		# concatenate g1 and g2 coords and convert to csv, and make header
		#chom1,start,stop,chrom2,start,stop
		echo '#chom1,start1,stop1,chrom2,start2,stop2' > ${anch}.g1.g2.coords
		paste -d , ${anch}.g1.coords ${anch}.g2.coords | sed -e 's/\s\+/,/g' >> ${anch}.g1.g2.coords

		# Writing .layout files (this section is not yet tested)
		layout=${x}.${y}.layout
		echo "Writing .layout files" && echo ""

		echo '# y, xstart, xend, rotation, color, label, va, bed' > ${layout}
		echo " .6,     .1,    .8,       0,      , ${x}, top, ${x}.bed" >> ${layout}
		echo " .4,     .1,    .8,       0,      , ${y}, bottom, ${y}.bed" >> ${layout}
		echo '# edges' >> ${layout}
		echo "e, 0, 1, ${x}.${y}.anchors.simple" >> ${layout}

		# Write .seqid files (this section is not yet tested)
		seqids=${x}.${y}.seqids
		rm -f ${seqids}
		echo "Writing .seqids files" && echo ""

		num=$(cat ${x}.bed | awk '{print $1}' | sort | uniq | wc -l | head -1)
		cat ${x}.bed | awk '{print $1}' | sort | uniq | awk '{printf $1 (!(NR % cols) ? c="\n" : c=" "); } END {print "";}' cols=${num} | sed 's/[[:blank:]]/,/g' | sed '/^$/d'> ${seqids}
		num=$(cat ${y}.bed | awk '{print $1}' | sort | uniq | wc -l | head -1)
		cat ${y}.bed | awk '{print $1}' | sort | uniq | awk '{printf $1 (!(NR % cols) ? c="\n" : c=" "); } END {print "";}' cols=${num} | sed 's/[[:blank:]]/,/g' | sed '/^$/d' >> ${seqids}

		# draw karyotype
		echo "Drawing karyotype.pdf" && echo ""
		python -m jcvi.graphics.karyotype ${seqids} ${layout}
		mv karyotype.pdf ${x}.${y}.karyotype.pdf
	done
done

# Remove self-self matches
for x in ${list[@]}; do
	rm -f ${x}.${x}.*
done
