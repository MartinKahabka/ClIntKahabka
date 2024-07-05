README file für die promoter variant analysis (Folder PROMOTER)

Goal:
Das Ziel des Scriptes ist es, mehrere VCF von unterschiedlichen Patienten mit gegebenen Promoterregionen zu vergleichen. Also ob die Patienten Varianten auf den gegebenen Promoterregionen besitzen

Generelle Info Script(Für Begriffe siehe unten):
Die variantPromoterRegion.py file analysiert jeweils eine nur Input-VCF (siehe Begriff)
Das Script wird über die runVariantProm.sh file ausgeführt. Dort wird der Ort definiert, in dem die alle zu analysierenden Input-VCF liegen.

Anwendung runVariantProm.sh
6 Parameter in folgender Reihenfolge sind anzugeben:

1. name: Identifyer dieses Durchlaufes (siehe Begriff). Unterschiedliche Ergebnisse werden so voneinander seperat getrennt gespeichert. Für einen Durchlauf ist dieser Identifyer für alle Outputfiles gleich
2. output_path: Form von "$output_path/name/..."
Definiert den Ordner indem die Outputfiles gespeichert werden. Falls der Ordner nicht vorhanden ist wir einer erstellt. 
3. input_path: Path zu dem Ordner in dem die Input-VCF liegen. Dabei werden alle Dateien gelesen die dem Pattern "*.vcf" folgen. Vorsicht! Es werden nicht nur Files in dem angegebenen Ordner, sondern auf in Unterordner beachtet.
4. promoter_path: Der Path zu der Datei indem die Informationen zu den Promoter liegen.
5. start: Definiert die Länge des Promoters, downstream von der TSS des Promoters. Ist der Wert also bei "1000", so werden 1000 Basen downstream vom Promoter als Teil der Promoterregion erachtet und dort wird nach Varianten gesucht.
6. end: Identisch mit 'start', allerdings upstream vom Promoter

Beispiel (von /kahabka/ClIntKahabka/promoter aus)
sh ./runVariantProm.sh "standart" "../../promoterResults" "../../../old/ClInt/combined_folder_dna/vcf_dna" "../../dEGenes/GRCh37_promoterChrPos.bed" "500" "100"  


Begriffe
Input-VCF: Die VC-file eines Patient
Durchlauf: Eine Ausführung des 'runVariantProm.sh' Scriptes
