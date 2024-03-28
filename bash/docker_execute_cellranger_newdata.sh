GRPHOME=~/path/to/grphome
PROJHOME=${GRPHOME}/path/to/prjhome
WORKHOME=/path/to/processed_data_new
IDS=$(id -u):$(id -g)
docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day3rep1 --csv=/prjhome/metadata/cellrangerconfig_G14_day3rep1.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day6rep1 --csv=/prjhome/metadata/cellrangerconfig_G15_day6rep1.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day10rep1 --csv=/prjhome/metadata/cellrangerconfig_G16_day10rep1.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day12rep1 --csv=/prjhome/metadata/cellrangerconfig_G17_day12rep1.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day14rep1 --csv=/prjhome/metadata/cellrangerconfig_G18_day14rep1.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day19rep1 --csv=/prjhome/metadata/cellrangerconfig_G19_day19rep1.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day28rep1 --csv=/prjhome/metadata/cellrangerconfig_G20_day28rep1.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day3rep2 --csv=/prjhome/metadata/cellrangerconfig_G21_day3rep2.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day6rep2 --csv=/prjhome/metadata/cellrangerconfig_G22_day6rep2.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day10rep2 --csv=/prjhome/metadata/cellrangerconfig_G23_day10rep2.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day12rep2 --csv=/prjhome/metadata/cellrangerconfig_G24_day12rep2.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day14rep2 --csv=/prjhome/metadata/cellrangerconfig_G25_day14rep2.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day19rep2 --csv=/prjhome/metadata/cellrangerconfig_G26_day19rep2.csv --localmem=90 --localcores=12 

docker run --name "${USER}_cellranger3" -it --rm --user ${IDS} -v ${WORKHOME}:/work -v ${GRPHOME}:/grphome -v ${PROJHOME}:/prjhome cancerbits/cellranger-7.1.0 multi --id=day28rep2 --csv=/prjhome/metadata/cellrangerconfig_G27_day28rep2.csv --localmem=90 --localcores=12 

