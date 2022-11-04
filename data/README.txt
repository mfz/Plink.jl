Test files

Started with a .fam, .bim, .bed triple.
Converted to plink.map, plink.ped using julia function Plink.make_ped

plink.fam, plink.bim, and plink.bed were created using the following command

plink1.9 --file plink --make-bed --out plink