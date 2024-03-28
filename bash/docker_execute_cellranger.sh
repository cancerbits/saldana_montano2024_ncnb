GRPHOME=/path/to/grphome
PROJHOME=${GRPHOME}/path/to/projhome
WORKHOME=/path/to/processed_data
IDS=$(id -u):$(id -g)
#rm -rf ${WORKHOME}  #uncomment if restarting from scratch
#mkdir -p ${WORKHOME} #uncomment if first time

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G1_GEX --csv=/prjhome/metadata/cellrangerconfig_G1_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G2_GEX --csv=/prjhome/metadata/cellrangerconfig_G2_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G3_GEX --csv=/prjhome/metadata/cellrangerconfig_G3_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G4_GEX --csv=/prjhome/metadata/cellrangerconfig_G4_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G5_GEX --csv=/prjhome/metadata/cellrangerconfig_G5_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G6_GEX --csv=/prjhome/metadata/cellrangerconfig_G6_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G7_GEX --csv=/prjhome/metadata/cellrangerconfig_G7_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G8_GEX --csv=/prjhome/metadata/cellrangerconfig_G8_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G9_GEX --csv=/prjhome/metadata/cellrangerconfig_G9_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G10_GEX --csv=/prjhome/metadata/cellrangerconfig_G10_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G11_GEX --csv=/prjhome/metadata/cellrangerconfig_G11_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G12_GEX --csv=/prjhome/metadata/cellrangerconfig_G12_GEX.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger2" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=G13_GEX --csv=/prjhome/metadata/cellrangerconfig_G13_GEX.csv --localmem=90 --localcores=12 

