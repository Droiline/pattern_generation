IF pbimx = 0 THEN print "pbimx for huge array is not defined":END
ilm = pbilm: imx = pbimx: imy = pbimy:
DIM huge ax(ilm, imx, imy)
PRINT "array huge "; imx
