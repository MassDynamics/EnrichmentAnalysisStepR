SPECIES?="human"
MSIG_VERSION="7.5"
BUILD_PACKAGE=False

install:
	Rscript install.R

go:
	./inst/scripts/run_go.sh ${SPECIES}

reactome:
	./inst/scripts/run_reactome.sh ${SPECIES}

msigdb:
	./inst/scripts/run_misgdb.sh ${MSIG_VERSION}

human:
	./inst/scripts/create_annotations.sh ${SPECIES} ${MSIG_VERSION}

all:
	./inst/scripts/make_all_annotations.sh ${BUILD_PACKAGE} ${MSIG_VERSION}
