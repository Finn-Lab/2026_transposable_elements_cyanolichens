set -ex 

SPECIES=(
	Leptogium_burgessi
	Lobaria_pulmonaria
	Peltigera_hymenina
	Peltigera_praetexata
	Ricasolia_virens
)

SHORT=(
	glLepBurg3
	glLobPulm2
	glPelHyme1
	glPelPrae3
	glRicVire11
)

for i in `seq 0 5`; do

	mkdir -p stringtie_combo_final/${SPECIES[$i]}
	
	awk '$7=="+"' stringtie_combo/${SPECIES[$i]}/${SPECIES[$i]}.fwdstringtie_v2.gtf > \
		stringtie_combo/${SPECIES[$i]}/${SPECIES[$i]}.fwdstringtie_v2.fwd_only.gtf
	
	awk '$7=="-"' stringtie_combo/${SPECIES[$i]}/${SPECIES[$i]}.revstringtie_v2.gtf > \
		stringtie_combo/${SPECIES[$i]}/${SPECIES[$i]}.revstringtie_v2.rev_only.gtf

	gffcompare -o stringtie_combo_final/${SPECIES[$i]}/${SHORT[$i]} \
		--strict-match -e 250 -d 100 -p ${SHORT[$i]} \
		stringtie_combo/${SPECIES[$i]}/${SPECIES[$i]}.fwdstringtie_v2.fwd_only.gtf \
		stringtie_combo/${SPECIES[$i]}/${SPECIES[$i]}.revstringtie_v2.rev_only.gtf

	gffread -g refs/combined_fasta/${SHORT[$i]}.fa \
		-w stringtie_combo_final/${SPECIES[$i]}/${SHORT[$i]}.tx.fa  stringtie_combo_final/${SPECIES[$i]}/${SHORT[$i]}.combined.gtf 
done

	
