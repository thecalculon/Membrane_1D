fol="example/data/"
pickled_data(fol)
xed,meany,meanperr,meanmerr=mean_indexp(fol)
xed,meany,meanperr,meanmerr=chop(xed,meany,meanperr,meanmerr,min_=0)