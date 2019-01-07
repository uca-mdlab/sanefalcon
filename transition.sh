#!/usr/bin/env bash

## give the correct directory as arg
ls -1 | while read line; do NAME=$(echo $line | cut -d '.' -f1-3,7) ; mv $line $NAME ;done

./combProf.sh /Users/davidpratella/projet_ff/sanefalcon/results_data_test_on_set_a_nucltrack_straver >> /Users/davidpratella/projet_ff/sanefalcon/nucleosome_profil/nucleosome_profile_set_a_straver/upstreamProfs1_set_a_straver.csv

./combProfI.sh /Users/davidpratella/projet_ff/sanefalcon/results_data_test_on_set_a_nucltrack_straver >> /Users/davidpratella/projet_ff/sanefalcon/nucleosome_profil/nucleosome_profile_set_a_straver/downstreamProfs1_set_a_straver.csv

## supprime la première colonne et reverse toutes les lignes (pour le fichier upstream)
awk -F "," '{$1=""; print $0}' upstreamProfs1_set_a_straver.csv | while read -ra words; do for ((i=${#words[@]}-1; i>=0; i--)); do printf "%s " "${words[i]}"; done; echo; done > upstreamProfs1_set_a_straver_bis_reverse.csv

## supprime la 1ere colonne (pour le fichier downstream)
awk -F "," '{$1=""; print $0}' downstreamProfs1_set_a_straver.csv > downstreamProfs1_set_a_straver_bis_reverse.csv

## paste la colonne avec les noms des samples puis avec le fichier upstream reversed et le fichier downstream
awk -F "," '{print $1}' downstreamProfs1_set_a_straver.csv | paste /dev/stdin upstreamProfs1_set_a_straver_bis_reverse.csv downstreamProfs1_set_a_straver_bis_reverse.csv > profile_nucleosome1_set_a_straver.txt