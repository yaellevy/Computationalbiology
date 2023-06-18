------------------
main.py

BS168.gb
main.py
Uniprot.xlsx
requirements.txt
י
install requirements.txt

--------------------------------------------------------------------------------------
Gene=Genes(11,'BS168.gb') – הפקודה הזו כבר קיימת בקובץ ולא מכוכבת – אין לככב!

                                                                                                                                                    -------------קבצים שישמרו במחשב לאחר ההרצה  ': parta.csv----------------------
  1Part ------------------------------------------------------------------------------
1.print(Gene.sortGene())
2.1 print(geneLength(Gene.features_list,'gene'))
2.2
print(geneLength(Gene.list_of_gene_cds,'gene'))
print(geneLength(Gene.list_of_gene_rna,'gene'))
print(gene_cds_statistics(Gene.list_of_gene_cds,'gene'))
2.3
print(gene_cds_statistics(Gene.list_of_gene_no_cds,'gene'))
print(Gene.gene_cds_statistics(Gene.list_of_gene_cds,'gene'))
2.4
histogram(geneLength(Gene.features_list, 'gene'), 'Length', 'Frequency', 200)
histogram(geneLength(Gene.list_of_gene_cds, 'gene'), 'Length', 'Frequency', 200)
histogram(geneLength(Gene.list_of_gene_no_cds, 'gene'), 'Length', 'Frequency', 200)

---------------------------------------------------------------------------------------------------------
3.1
print('mean of %AT in all Genom:',countAT(Gene.sequence))#mean %AT in all Genom
3.2
print('mean of %AT in Genes that coding to CDS:',
     getMeanValue(Gene.countGenesAT(Gene.list_of_gene_cds, 'gene')))  # mean of %AT in Genes that coding to CDS
print('mean of %AT in Genes that Are Not Coding to CDS:', getMeanValue(
    Gene.countGenesAT(Gene.list_of_gene_no_cds, 'gene')))  # mean of %AT in Genes that Not coding to CDS

3.4
histogram(Gene.countGenesAT(Gene.list_of_gene_cds,'gene'),'Length','Frequency',200)
3.5
print(Gene.getDictFiveTopFiveLess())
---------------------------------------------------------------------------------------------------------
4.1
count, llistGen = Gene.cellWall()
cellWallDict = geneLength(llistGen, 'CDS')
print('num of cellWal
4.2
count,llistGen = Gene.cellWall()
cellWallDict   = geneLength(llistGen,'CDS')
print(cellWallDict)
print(gene_cds_statistics(llistGen,'CDS'))
histogram(cellWallDict,'Length Cell Wall Genes','Frequency',200)
4.3
count, llistGen = Gene.cellWall()
print(Gene.countGenesAT(llistGen, 'CDS'))
histogram(Gene.countGenesAT(llistGen, 'CDS'), '%AT For Cell Wall Genes', 'Frequency', 200)
print(At_statistics(Gene.countGenesAT(llistGen, 'CDS')))
---------------------------------------------------------------------------------------------------------
5
                                                                                                                                ------------קבצים שישמרו במחשב לאחר הרצה  ': gene_exceptions.csv-------------------------------              
Gene.compareGenesCds()



2.1
Gene.uniProtCompare()
2.2.1
dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
print(transmemranStatistics(dictMamber))
histogram(dictMamber,'Length','Frequency',10)	
2.2.2
dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
print(transmemranStatistics(aminoacid))
histogram(aminoacid,'Length','Frequency',10)
2.3.1
dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
histogram(dicB,'%AT','Frequency',10)
2.3.2
dictMamber,dictSeq,aminoacid,dicA,dicB,dicC = Gene.getTransmambernal()
ATPercenrHistogram(dicA,dicB,dicC)
print(At_statistics(dicB))
print(At_statistics(dicC))



                                                                                                                                                                      שני קבצי GB בשם: ["CovJune21.gb", "CovDecember22.gb"]    


רעיון החלוקה של חלק ג' כל הפונקציות עם המשתנים הרלוונטים מצויות בראש כל חלק,
ומתחתיו שורות קוד הרצה למענה על השאלות.
ניתן להריץ את כל הקובץ על ידי run וכל השורות יודפסו עם התשובות אחת אחרי השנייה.
שורות ההרצה מפורטות להלן:

###------------SECTION 1 Synonymous Postions Per Codon------------###
הרצת החלק הראשון שעונה על שאלה 1 ניתן להריץ בתחתית בלוק הפונקציות של חלק זה 
בשורת קוד מספר: 66 יודפס המילון עם התוצאות.
print("SECTION 1 Answer:", '\n', count_synonymous_positions_for_codon())

במידע ויש רצון לייצא כקובץ csv יש לבטל את הערות בשורת קוד מתחת ולייצא.
 df = pd.DataFrame(count_synonymous_positions_for_codon().items())
df.to_csv('section1_csv.csv', index=False, header=False)

###------------SECTION 2 Comparing The Two Viruses------------###
הרצת section 2 מתחילה בהזנת מייל וקריאת קבצי הGB שיירדו בעזרת הפונקציות החל משורה 213

המענה על שאלה 2. א ניתן למצוא תחת הכותרת
# --Section 2.a: Answer for how many genes and proteins there are in the viruses--#
על ידי הרצת שורות קוד 229,231:
print_numbers_of_genes_proteins(accessions[0], count_gene_and_cds(record_cov_june21)[0],
                                count_gene_and_cds(record_cov_june21)[1])
print_numbers_of_genes_proteins(accessions[1], count_gene_and_cds(record_cov_december22)[0],
                                count_gene_and_cds(record_cov_june21)[1])

המענה על שאלה 2. ב ניתן למצוא תחת הכותרת
# --Section 2.b: Answer for how many genes are common between the viruses--#
על ידי הרצת שורת קוד 236:
print(common_genes(record_cov_june21, record_cov_december22))

המענה על שאלה 2. ג ניתן למצוא תחת הכותרת
# --Section 2.c: Answer for dnds for 5 chosen genes--#
מבוצע תהליך בחירת הגנים ובדיקת תרגומם החל משורה 240:
# chosen list of genes
chosen_genes_compare = ['ORF3a', 'E', 'M', 'N', 'ORF6']

בשורה 253 ישנה ללואה שמריצה את הפונקציות הרלוונטיות לחישוב מדד dnds
ומדפיס את התוצאות עבור החלבונים הרצויים.
print("Section 2.c: Answer:", '\n')
for gene in chosen_genes_compare:

