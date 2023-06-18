from main import  Gene


def test_get_transmambernal():
   assert (('gene' in Gene.sortGene().keys()) == True) and (('CDS' in Gene.sortGene().keys()) == True)

