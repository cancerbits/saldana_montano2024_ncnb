vl=fread(file="~/ncnb2_code/metadata/velocyto_libraries.csv")

lib2G= Vectorize(function(x) dayrep.dataset.map[[x]], USE.NAMES=F)

vl2=vl %>% mutate(orig.ident= lib2G(library_name))

fwrite(vl2, file="~/ncnb2_code/metadata/velocyto_libraries_withident.csv")