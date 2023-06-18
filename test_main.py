from main import  Gene


def test_get_transmambernal():
   assert ('tRNA' in Gene.sortGene().keys()) == True
